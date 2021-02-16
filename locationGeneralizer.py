########################################################################
#  Copyright 2020 Battelle Energy Alliance, LLC   ALL RIGHTS RESERVED  #
#  Mobility Systems & Analytics Group, Idaho National Laboratory       #
########################################################################

# Location Generalizer
# Release 1.0 2/16/2021

import pyodbc
import pandas as pd
import pickle
from datetime import datetime, timedelta
import time
import math
import yaml
from pathlib import Path
import csv
import numpy as np
from sklearn.cluster import DBSCAN
from shapely import geometry
from shapely.geometry import MultiPoint
from haversine import haversine, Unit
import pynput
from pandasql import sqldf

class cfg():
    with open('locationGeneralizer.yml') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    
    odbcConnectionString=config['odbcConnectionString']
    inputTableOrCSV= config['inputTableOrCSV']
    vehiclesInChunk = config['vehiclesInChunk']

    qryVehicleIDList =config['qryVehicleIDList']
    qryVehicleInfo = config['qryVehicleInfo']
    qryVehicleIDList = qryVehicleIDList.replace('{inputsrc}', inputTableOrCSV)
    qryVehicleInfo = qryVehicleInfo.replace('{inputsrc}', inputTableOrCSV)

    errorLogFileName = config['errorLogFileName']
    heartbeatFileName = config['heartbeatFileName']
    locationInfoFileName = config['locationInfoFileName']
    homeInfoFileName = config['homeInfoFileName']
    pklCensusDivisionsFileName = config['pklCensusDivisionsFileName']
    evseLookupFileName = config['evseLookupFileName']

    bboxes = config['boundingBoxes']

    gpsOdoThreshold_mi = config['gpsOdoThreshold_mi']
    minTrips = config['minTrips']
    minLastTrips = config['minLastTrips']
    minPctParks = config['minPctParks']
    numL2Rounding = config['numL2Rounding']
    numDCRounding = config['numDCRounding']

    doCheck = config['doCheck']

    dayEndHours = config['dayEndHours']
    dayEndMinutes = config['dayEndMinutes']

    dbscan_eps_ft = config['dbscan_eps_ft']
    dbscan_min_spls = config['dbscan_min_spls']

    evseDistRange_Miles = config['evseDistRange_Miles']
    evseLatRange = config['evseLatRange']
    evseLonRange = config['evseLonRange']

    hdrErrorLogCSV = config['hdrErrorLogCSV']
    hdrLocationInfoCSV = config['hdrLocationInfoCSV']
    hdrHomeInfoCSV = config['hdrHomeInfoCSV']
    colLocationInfo = config['colLocationInfo']
    colHomeInfo = config['colHomeInfo']

    verbosity = config['verbosity']
    stopProcessing = False

    errFilePath = Path(errorLogFileName)
    if not errFilePath.exists():
        # ErroLog output file
        hdr = pd.DataFrame(hdrErrorLogCSV)
        hdr.to_csv(errorLogFileName, index=False, header=False, mode='w')
    # use one line buffering - every line written is flushed to disk
    errorFile = open(errorLogFileName, mode='a', buffering=1, newline='')
    errorWriter = csv.writer(errorFile)

def main():
    # trust chained assignments (no warnings)
    pd.set_option('mode.chained_assignment', None)

    # LocationInfo output file
    locationFilePath = Path(cfg.locationInfoFileName)
    if not locationFilePath.exists():
        hdr = pd.DataFrame(cfg.hdrLocationInfoCSV)
        hdr.to_csv(cfg.locationInfoFileName, index=False, header=False, mode='w')

    # HomeInfo output file
    homeFilePath = Path(cfg.homeInfoFileName)
    if not homeFilePath.exists():
        hdr = pd.DataFrame(cfg.hdrHomeInfoCSV)
        hdr.to_csv(cfg.homeInfoFileName, index=False, header=False, mode='w')

    ## geopandas can read the shapefile directly, but we pickled it into one file
    ## a single pickle file simplifies distribution whereas,
    ## loading a shapefile requires several adjacent accompanying files
    divisions = pd.read_pickle(cfg.pklCensusDivisionsFileName)

    # get Public EVSE stations
    EVSEs = pd.read_csv(cfg.evseLookupFileName)

    # pyodbc attempts to turn off autocommit before returning connection and this causes file connections (like CSVs), which do not support transactions, to fail
    # when autocommit is explictly set, pyodbc will not attempt any changes
    cnxn = pyodbc.connect(cfg.odbcConnectionString, autocommit=True)

    lastVehicle = 0
    hbFilePath = Path(cfg.heartbeatFileName)
    if hbFilePath.exists():
        with open(hbFilePath, 'r') as hb:
            lastVehicle = hb.readline()
        cfg.errorWriter.writerow([datetime.now(), lastVehicle, -1,'Restarting after vehicle {}'.format(lastVehicle)])
        print('Restarting after vehicle {}'.format(lastVehicle))

    # get sorted list of all vehicle IDs
    qry = cfg.qryVehicleIDList.replace('{startVehicle}', str(lastVehicle))
    df = pd.read_sql(qry, cnxn)
    numOfVehicles = cfg.vehiclesInChunk # number of vehicle to process at a time. We can't process all at once due to dataset size, so this is the "chunk size" to process
    vehicleList =  df['VehicleID'].tolist()
    # divide up vehicle ID list into sections of <numOfVehicle> length chunks (we'll read data in one chunk at a time to avoid memory overrun)
    chunks = [vehicleList[i * numOfVehicles:(i+1)*numOfVehicles] for i in range((len(vehicleList) + numOfVehicles -1) // numOfVehicles)]
    i = 0
    vcnt = 0
    for chunk in chunks:
        chunkList = ','.join(str(e) for e in chunk)
        qry = cfg.qryVehicleInfo.format(chunkList) # insert vehicleIDs into "in" list
        if cfg.verbosity > 0: print('Fetching data')
        chunkData = pd.read_sql(qry, cnxn, parse_dates=['TripStartLocalTime', 'TripEndLocalTime'])
        # create new column for flag to exclude bad records
        chunkData['Include'] = True
        i += 1
        print("chunk: {}, vehicle from {} through {}".format(i, chunk[0], chunk[-1]))

        # iterate through one vehicle at a time
        for v in chunk:
            if cfg.stopProcessing: exit()
            if cfg.verbosity > 0: print('Vehicle: {}'.format(v))

            # create empty LocationInfo data frame
            # GPS coordinates are added here for convenience, but will not be carried into LocationInfo output file
            locationInfo = pd.DataFrame(columns = cfg.colLocationInfo)
            
            # create empty HomeInfo data frame
            homeInfo = pd.DataFrame(columns = cfg.colHomeInfo)
            homeInfo['HomeStartLocalTime'] = pd.NaT
            homeInfo['HomeEndLocalTime'] = pd.NaT            
            vcnt += 1
            
            # grab all records in vehicle v
            vData = chunkData[chunkData['VehicleID'] == v]
            # create new column to check for Odometer gaps, i.e missing trips
            vData['resid_Miles'] = vData['TripStartOdometer_Miles'].shift(periods=-1) - vData['TripEndOdometer_Miles']

            ### Check validity of data, marking invalid records (Include = True/False)
            if cfg.verbosity > 1: print(' Check for valid values')
            if cfg.doCheck: vData = DoValidityChecking(v, vData)
            vData.resid_Miles = vData.resid_Miles.astype(object).where(vData.resid_Miles.notnull(), None) # set NaN to None (becomes Null for DB)
            # toss out rows that failed vailidity check
            vData = vData[vData.Include == True]
            numTrips = len(vData)
            if numTrips < cfg.minTrips:
                if cfg.verbosity > 1: print(' Not enough trips, vehicle skipped.')
                cfg.errorWriter.writerow([datetime.now(), v, -1,'Not enough trips, vehicle skipped. ({} need >= {})'.format(numTrips, cfg.minTrips)])
            else:
                # create new column for identify first/last trip of day
                vData['TripFlag'] = None
                ### Identify first and last of trip of day
                if cfg.verbosity > 1: print(' Defining first/last trip of day')
                vData = flagTrips(v, vData)
                ### Find clusters of vehicle locations
                if cfg.verbosity > 1: print(' Clustering')
                vData = clusterData(v, vData)

                # # drop rows - remove previous vehicle info
                # homeInfo.drop(homeInfo.index, inplace=True)
                # locationInfo.drop(locationInfo.index, inplace=True)

                # add row to LocationInfo data frame
                liList = [vData[['VehicleID', 'TripStartLocalTime', 'TripEndLocalTime', 'TripStartLatitude', 'TripStartLongitude', 'TripEndLatitude','TripEndLongitude', 'TripStartClusterID', 'TripEndClusterID']]]
                locationInfo = locationInfo.append(liList, ignore_index=True)
            
                # add ParkEndLocalTime for convenience - its the same as TripStartLocalTime in next row
                vData['ParkEndLocalTime'] = vData['TripStartLocalTime'].shift(periods=-1)
                vData['ParkDuration_hr'] = (vData['ParkEndLocalTime'] - vData['TripEndLocalTime'])/np.timedelta64(1,'h')

                ########################
                #### FIND HOME DECISION TREE return more than one home, but must returns an array to conform with other methods that may return more that one)
                if cfg.verbosity > 1: print(' Identifying home location')
                vData = findHome_DecisionTree(v, vData)
                homeClusters = list(set(vData[vData['location'] == 'home']['TripEndClusterID']))
                ########################

                if len(homeClusters) == 0:
                    cfg.errorWriter.writerow([datetime.now(), v, -1,'No home clusters found - vehicle skipped.'])
                    continue # continue with next vehicle

                ########################
                #### PROCESS HOME AND LOCATION INFO returns the data we want to push as output files
                if cfg.verbosity > 1: print(' Calculating output metrics')
                isOK, locationInfo, homeInfo = processHome(v, divisions, vData, locationInfo, homeInfo, homeClusters, EVSEs)
                if not isOK:
                    continue
                ########################

                ############################# IMPORTANT #########################################
                #  CLEANUP HONEINFO AND LOCATIONINFO FOR EXPORT (remove sensitive data)
                homeInfo.drop(homeInfo[homeInfo['Primary'].isnull()].index, inplace=True)
                homeInfo.drop(['CentroidLatitude', 'CentroidLongitude', 'Primary'], axis=1, inplace=True)
                locationInfo.drop(['TripStartLatitude','TripStartLongitude','TripEndLatitude','TripEndLongitude','TripStartClusterID','TripEndClusterID'], axis=1, inplace=True)
                #################################################################################

                # write to output files
                if cfg.verbosity > 1: print(' Writing to output files')
                locationInfo.to_csv(cfg.locationInfoFileName, index=False, header=False, mode='a')
                homeInfo.to_csv(cfg.homeInfoFileName, index=False, header=False, mode='a')

            # # use one line buffering - every line written is flushed to disk
            with open(cfg.heartbeatFileName, mode='w', buffering=1, newline='') as hb:
                hb.write(str(v))

def findHome_DecisionTree(v, vData):
    vData['location'] = 'unknown'

    ## apply filter rules to get qualifying clusters
    vQual = rule1(vData) # this does not find good clusters, but checks if there is enough vehicle data to continue
    if vQual.empty: 
        cfg.errorWriter.writerow([datetime.now(), v, -1,'Vehicle failed at rule 1'])
        return vData # leave all clusters marked unknown
    
    vQual = rule2(vData) # vData (the complete datset) is passed and we start eliminating clusters to find only those of interest
    if vQual.empty:  
        cfg.errorWriter.writerow([datetime.now(), v, -1,'Vehicle failed at rule 2'])
        return vData # leave all clusters marked unknown
    
    vQual = rule3(vQual) # vQual is passed and further filtered
    if vQual.empty:
        cfg.errorWriter.writerow([datetime.now(), v, -1,'Vehicle failed at rule 3'])
        return vData # leave all clusters marked unknown

    qualClusters = vQual['TripEndClusterID']

    vData['location'] = 'away'
    # apply decision tree to get home clusters
    homeClusters = decision1(vData, qualClusters) # find home by 65% rule
    if not homeClusters.empty:
        vData.loc[vData.TripEndClusterID.isin(homeClusters.TripEndClusterID), ['location']] = 'home'

    homeClusters = decision2(vData, qualClusters) # look for home in remaining clusters
    if not homeClusters.empty:
        vData.loc[vData.TripEndClusterID.isin(homeClusters.tolist()), ['location']] = 'home'

    return vData

def rule1(vQual):
    # 1. vehicle must have at least 30 days of known last trips
    ids = vQual[((vQual['TripFlag'] == 'L') | (vQual['TripFlag'] == 'FL')) & (vQual['resid_Miles'] > -1) & (vQual['resid_Miles'] < 1)]['TripEndClusterID']
    if len(ids) < 30:
        vQual = vQual[0:0]
    return vQual

def rule2(vQual):
    # 2. get clusters with > 21 days between first and last parks (i.e. 'cluster period')
    #### get number of days, grouped by cluster id, between first and last park
    ids = vQual.groupby('TripEndClusterID')['TripEndDateOffset'].max() - vQual.groupby('TripEndClusterID')['TripEndDateOffset'].min()
    if ids.empty:
        vQual = vQual[0:0]
    else:
        ids = ids.dt.days
        #### add one day to max-min range
        ids = ids.add(1)
        #### get cluster ids with more than 21 days between first and last park
        ids = ids[ids > 21]
        if ids.empty:
            vQual = vQual[0:0]
        else:
            #### get data belonging to qualifying clusters
            vQual = vQual[vQual['TripEndClusterID'].isin(ids.index.tolist())]
    return vQual

def rule3(vQual):
    # 3. get clusters with >= 15 parks during the cluster period
    #### get number of parks, grouped by cluster id
    ids = vQual.groupby('TripEndClusterID')['TripEndClusterID'].count()
    #### get cluster ids having at least 15 parks
    ids = ids[ids >= 15]
    if ids.empty:
        vQual = vQual[0:0]
    else:
        #### get data belonging to qualifying clusters
        vQual = vQual[vQual['TripEndClusterID'].isin(ids.index.tolist())]
    return vQual

def decision1(vData, qualClusters):
    # how many days in period, grouped by cluster
    ids = vData[((vData['TripFlag'] == 'L') | (vData['TripFlag'] == 'FL')) & (vData['resid_Miles'] > -1) & (vData['resid_Miles'] < 1)][['TripEndClusterID','TripEndDateOffset']]
    daysInPeriod = ids.groupby('TripEndClusterID').count()  # NumDaysInPeriodWithParkAKLTInCluster 
    
    # how many driving days in period, grouped by cluster
    vMin = vData.groupby('TripEndClusterID')['TripEndDateOffset'].min()
    vMax = vData.groupby('TripEndClusterID')['TripEndDateOffset'].max()
    vRange = pd.concat([vMin, vMax], axis=1)
    vRange.columns=['minDate', 'maxDate']

    q = '''
    select TripEndClusterId, count(*) numDrivingDays 
    from (
        select distinct v.TripEndClusterID, v.minDate, i.TripEndDateOffset, v.maxDate 
        from ids i 
        inner join vRange v
        on i.TripEndDateOffset >= v.minDate 
        and i.TripEndDateOffset <= v.maxDate
    ) a 
    group by TripEndClusterID
    '''
    drivingDays = sqldf(q, locals()) # NumDrivingDaysInPeriodWithKnownLastTrip

    # percent days with park, grouped by cluster
    q = '''
    select a.TripEndClusterID, a.TripEndDateOffset*1.0 / b.numDrivingDays parkDays_pct
    from daysInPeriod a
    inner join drivingDays b
    on a.TripEndClusterID = b.TripEndClusterID
    '''
    parkDaysPct = sqldf(q, locals()) # PercDaysWithKLTWithParkAKLT

    homeClusters = parkDaysPct[parkDaysPct['parkDays_pct'] > 0.65]
    homeClusters = homeClusters[homeClusters['TripEndClusterID'].isin(qualClusters.to_list())]
    return homeClusters

def decision2(vData, qualClusters):
    parkDistinctTimes = set(vData['TripEndDateOffset'])  #### matches check query

    missingTrips = set(vData[(vData['resid_Miles'] > 1) | (vData['resid_Miles'] < -1)]['TripEndDateOffset']) #### matches check query

    goodTimes = parkDistinctTimes.difference(list(missingTrips)) #### matches check query
    goodTimes = pd.DataFrame(list(goodTimes), columns=['TripEndDateOffset']) # convert set to dataframe

    vDataGood = vData[vData['TripEndDateOffset'].isin(goodTimes['TripEndDateOffset'])]
    totalParkedTimeHr = vDataGood.groupby('TripEndClusterID')['ParkDuration_hr'].sum()

    vMin = vData.groupby('TripEndClusterID')['TripEndDateOffset'].min()
    vMax = vData.groupby('TripEndClusterID')['TripEndDateOffset'].max()
    vRange = pd.concat([vMin, vMax], axis=1)
    vRange.columns=['minDate', 'maxDate']

    qry = '''
    Select A.TripEndClusterID, COUNT(*) NumDrivingDaysWithNoMissingTripsInPeriod
    from 
    ( Select distinct B.TripEndClusterID, B.MinDate, A.TripEndDateOffset, B.MaxDate
    from goodTimes A
    inner join vRange B
    on A.TripEndDateOffset >= B.minDate
    and A.TripEndDateOffset <= B.maxDate
    ) A
    group by TripEndClusterID
    '''
    numDrivingDays  = sqldf(qry, locals()) # NumDrivingDaysWithNoMissingTripsInPeriod

    totalParkedTimeHr = pd.DataFrame(totalParkedTimeHr)
    qry = '''
    select b.TripEndClusterID, ParkDuration_hr /  NumDrivingDaysWithNoMissingTripsInPeriod TotalParkedTimeInPeriodInCluster_hr_PerDrivingDayWNMT 
    from totalParkedTimeHr a
    inner join numDrivingDays b
    on a.TripEndClusterID = b.TripEndClusterID
    '''
    totalParked = sqldf(qry, locals())
    homeClusters = totalParked[totalParked['TotalParkedTimeInPeriodInCluster_hr_PerDrivingDayWNMT'] > 9]['TripEndClusterID']
    homeClusters = homeClusters[homeClusters.isin(set(qualClusters))]

    return homeClusters

def getEVSEDistance(row, homeLat, homeLong):
    dist = haversine((row.Latitude, row.Longitude), (homeLat, homeLong), unit=Unit.MILES)
    return dist

def getStartLocationDistance(row, homeLat, homeLong, homeStart, homeEnd):
    if (homeStart <= row['TripEndLocalTime'] <= homeEnd):
        startDist = round(haversine((row['TripStartLatitude'], row['TripStartLongitude']), (homeLat, homeLong), unit=Unit.MILES))
    else:
        startDist = row['TripStartDistanceFromHome_Miles']
    return startDist

def getEndLocationDistance(row, homeLat, homeLong, homeStart, homeEnd):
    if (homeStart <= row['TripEndLocalTime'] <= homeEnd):
        endDist = round(haversine((row['TripEndLatitude'], row['TripEndLongitude']), (homeLat, homeLong), unit=Unit.MILES))
    else:
        endDist = row['TripEndDistanceFromHome_Miles']
    return endDist

def getLocationInfoData(row, HomeInfo, HomeClusterID):
    return 1

def isInHomeCluster(se, TripClusterID, TripLocaltTime, homeInfo):
    homeID = -1
    inRange = False
    ## is park in a home cluster (i.e. is park's cluster ID equal to a home cluster ID)
    if TripClusterID in set(homeInfo['HomeID']):
        homeID = TripClusterID
        ## in a home cluster, but is park time in home cluster period
        inRange = isInRangeSet(se, homeInfo[homeInfo['HomeID'] == homeID][['HomeStartLocalTime', 'HomeEndLocalTime']], TripLocaltTime)
            #inRange = True
    return homeID, inRange

def isInRangeSet(se, homeStartEnd, locationTime):
    if se == 'Start':
        if len(homeStartEnd[(homeStartEnd['HomeStartLocalTime'] < locationTime) & (homeStartEnd['HomeEndLocalTime'] >= locationTime)]) > 0:
            return True
        return False
    else:
        if len(homeStartEnd[(homeStartEnd['HomeStartLocalTime'] <= locationTime) & (homeStartEnd['HomeEndLocalTime'] > locationTime)]) > 0:
            return True
        return False

def isInRange(se, homeStartEnd, locationTime):
    if se == 'Start':
        if homeStartEnd['HomeStartLocalTime'] < locationTime <= homeStartEnd['HomeEndLocalTime']:
            return True
        return False
    else:
        if homeStartEnd['HomeStartLocalTime'] <= locationTime < homeStartEnd['HomeEndLocalTime']:
            return True
        return False

def isInTupleRange(se, start, end, locationTime):
    if se == 'Start':
        if start < locationTime <= end:
            return True
        return False
    else:
        if start <= locationTime < end:
            return True
        return False

def youHome(row, se, locationTime):
    if se == 'Start':
        if row['HomeStartLocalTime'] < locationTime <= row['HomeEndLocalTime']:
            row['locIn'] = True
        row['locIn'] = False
    else:
        if row['HomeStartLocalTime'] <= locationTime < row['HomeEndLocalTime']:
            row['locIn'] = True
        row['locIn'] = False


def processHome(v, divisions, vData, vLocationInfo, homeInfo, homeClusters, EVSEs):
    for cID in homeClusters:
        dfPts = vData[vData['TripEndClusterID'] == cID][['TripEndLatitude', 'TripEndLongitude']]
        mpPts = MultiPoint(dfPts.to_numpy())
        CP = mpPts.centroid
        CP = geometry.Point(CP.y, CP.x)
        for i, division in divisions.iterrows():
            if division.geometry.contains(CP):
                st = EVSEs[(EVSEs['Latitude'] > (CP.y - cfg.evseLatRange)) & 
                            (EVSEs['Latitude'] < (CP.y + cfg.evseLatRange)) & 
                            (EVSEs['Longitude'] > (CP.x - cfg.evseLonRange)) & 
                            (EVSEs['Longitude'] < (CP.x + cfg.evseLonRange))]
                if not st.empty:
                    st['hMiles'] = st.apply(getEVSEDistance, args=(CP.y, CP.x), axis=1)
                    st = st[st['hMiles'] <= cfg.evseDistRange_Miles]
                l2Cnt = 0
                dcCnt = 0
                if not st.empty:
                    l2Cnt = st['L2'].sum()
                    dcCnt = st['DCFC'].sum()
                    l2Cnt = round(l2Cnt, cfg.numL2Rounding)
                    dcCnt = round(dcCnt, cfg.numDCRounding)
                    if l2Cnt == 0: l2Cnt = 1
                    if dcCnt == 0: dcCnt = 1
                    homeStart = vData[vData['TripEndClusterID'] == cID]['TripEndLocalTime'].min()  
                    homeEnd = vData[vData['TripEndClusterID'] == cID]['TripEndLocalTime'].max()
                    
                    newRow = {'VehicleID':int(v), 'HomeID':cID,
                            'HomeStartLocalTime':homeStart, 'HomeEndLocalTime':homeEnd, 
                            'HomeRegion':division['NAME'], 'PublicChargingDensityL2':l2Cnt, 'PublicChargingDensityDCFC':dcCnt,
                            'CentroidLatitude':CP.y, 'CentroidLongitude':CP.x, 'Primary': False}
                    homeInfo = homeInfo.append(newRow, ignore_index=True)
                break # exit the division loop
        cfg.errorWriter.writerow([datetime.now(), v, -1,'No census division found for cluster.'])

    if homeInfo.empty:
        cfg.errorWriter.writerow([datetime.now(), v, -1,'No usable clusters found for homeInfo - vehicle skipped.'])
        return False, vLocationInfo, homeInfo


    # vehicles with only one home are marked as the primary home and primary home detection below is skipped
    if len(homeInfo) == 1:
        homeInfo['Primary'] = True
    else:
        # collect period start and period end date for each cluster
        # get cluster period start dates of each range
        sranges = vData[vData['TripEndClusterID'].isin(homeInfo['HomeID'])].groupby('TripEndClusterID')[['TripEndDateOffset', 'TripEndLocalTime']].min()
        eranges = vData[vData['TripEndClusterID'].isin(homeInfo['HomeID'])].groupby('TripEndClusterID')[['TripEndDateOffset', 'TripEndLocalTime']].max()
        sranges['period'] = 's'
        eranges['period'] = 'e'
        # assemble period start/end dates into sorted list of dates
        ranges = sranges.append(eranges)
        numDates = len(ranges['TripEndDateOffset'])
        numUniqDates = len(set(ranges['TripEndDateOffset']))
        if numDates != numUniqDates:
            cfg.errorWriter.writerow([datetime.now(), v, -1,'Range dates are not unique - vehicle skipped.'])
            return False, vLocationInfo, homeInfo

        ranges = ranges.sort_values(by=['TripEndDateOffset'])
        # make date ranges from first date to second date, second date to third date, etc.
        rangesEnd =  ranges.shift(-1)
        rangesEnd.rename(columns={'TripEndDateOffset': 'End', 'TripEndLocalTime': 'HomeEnd', 'period': 'endperiod'}, inplace=True)
        ranges.rename(columns={'TripEndDateOffset': 'Start', 'TripEndLocalTime': 'HomeStart', 'period': 'startperiod'}, inplace=True)
        ranges = pd.concat([ranges, rangesEnd], axis=1)
        ranges = ranges[:-1] # remove last row
        ranges = ranges.reset_index()

        # when a range start date originates from a cluster period end, it should be incremented
        eidxs = ranges[ranges['startperiod']=='e'].index
        erows = ranges[ranges.index.isin(eidxs)]
        erows['Start'] += timedelta(days=1)
        ranges[ranges.index.isin(eidxs)] = erows

        # when a range end date originates from a cluster period start, should be decremented
        sidxs = ranges[ranges['endperiod']=='s'].index
        srows = ranges[ranges.index.isin(sidxs)]
        srows['End'] += timedelta(days=-1)
        ranges[ranges.index.isin(sidxs)] = srows

        # create column for each home cluster for every park "after known last trip of day" (AKLT) days parked
        for x in list(homeInfo['HomeID']): ranges[x] = 0
        ranges.drop(['TripEndClusterID'], axis=1, inplace=True)

        # get cluster id and end date (offset) for every AKLT
        aklt = vData[((vData['TripFlag'] == 'L') | (vData['TripFlag'] == 'FL')) & (vData['resid_Miles'] > -1) & (vData['resid_Miles'] < 1)][['TripEndClusterID','TripEndDateOffset']]
        # loop through ranges, counting days parked in each cluster
        for ridx, rng in ranges.iterrows():
            st = rng['Start']
            en = rng['End']
            # get number of days parked in given range for each home cluster (returned as a Series)
            akltrow = aklt[(aklt['TripEndClusterID'].isin(homeInfo['HomeID'])) & (aklt['TripEndDateOffset'] >= st)  & (aklt['TripEndDateOffset'] <= en)].groupby('TripEndClusterID')['TripEndClusterID'].count()
            # write days parked to ranges dataframe
            for cid, numdays in akltrow.iteritems():
                rng[cid] = numdays
                ranges.iloc[ridx] = rng

        # initialize dataframes
        homeInRange = pd.DataFrame(columns = ['cID'])
        homeNumDays = pd.DataFrame(columns = ['cID'])
        for idx in ranges.index:
            homeInRange[idx] = 0
            homeNumDays[idx] = 0

        for cID in homeInfo['HomeID']:
            homeNumDays.loc[len(homeNumDays)] = [cID] + list(ranges[cID])

        # find ranges that are within each cluster's period
        for idx, homeID in homeInfo['HomeID'].iteritems():
            # initialize cluster's within range flag to 0
            homeInRange.loc[len(homeInRange.index)] = [homeID] + ([0] * len(ranges))
            # index list of ranges within homeID's period
            rin = list(ranges[(ranges['HomeStart'] >= homeInfo.iloc[idx]['HomeStartLocalTime']) & (ranges['HomeEnd'] <= homeInfo.iloc[idx]['HomeEndLocalTime'])].index)
            # set flag to 1 for ranges that are within the homeID's period
            row = homeInRange[homeInRange['cID'] == homeID]
            row[rin] = 1
            homeInRange[homeInRange['cID'] == homeID] = row

        # check number of cluster within range - no cluster = primary home, one cluster = cID is primary home, else leave for number of days check
        # get row indexes of homeInfoRange that have more then 1 home in range
        r = homeInRange.iloc[:,1:].sum(axis=0)
        multiInRangeidxs = r[r[0:] > 1]

        ranges['Primary'] = False
        # set ranges with single home to primary
        primaryIdxs = homeInRange.sum(axis=0)
        primaryIdxs = primaryIdxs[primaryIdxs==1]
        for idx in primaryIdxs.index:
            homeIDIdx = homeInRange[idx].idxmax(axis=0)
            homeID = homeInRange.iloc[homeIDIdx]['cID']
            row = ranges.iloc[idx]
            row['Primary'] = homeID
            ranges.iloc[idx] = row

        # in ranges with multi homes, if primary home can be determine, set it
        for col in list(multiInRangeidxs.index):
            # if range has multi homes set those with a single max num of days to primary, else not primary
            if len(homeNumDays[col][homeNumDays[col] == homeNumDays[col].max()]) == 1:
                idx = homeNumDays[col][homeNumDays[col] == homeNumDays[col].max()].index
                homeID = homeNumDays.iloc[idx]['cID'].item()
                row = ranges.iloc[col]
                row['Primary'] = homeID
                ranges.iloc[col] = row

        # create an add-on to homeInfo of ranges showing which HomeID was the primary home in that range
        newHomeInfo = pd.DataFrame()
        for i, row in ranges.iterrows():
            # copy the homeInfo row as starting point, then update the other fields in the new row
            nh = homeInfo[homeInfo['HomeID'] == row['Primary']]
            nh['HomeStartLocalTime'] = row['HomeStart']
            nh['HomeEndLocalTime'] = row['HomeEnd']
            nh['Primary'] = True
            newHomeInfo = newHomeInfo.append(nh, ignore_index=True)

        homeInfo['Primary'] = None
        homeInfo = homeInfo.append(newHomeInfo)
        homeInfo = homeInfo.reset_index(drop=True)

    ### Update vLocationInfo
    for se in ['Start', 'End']:
        # This mapping of column names to column index is nasty, but helps with performance.
        # Using the itertuples construct is much faster and is now used instead of iterrows. Iterrows
        # would allow for easire, more inteligible coding but tuples are more performant. Using
        # iterrows does not allow for so the column indexes must change instead.
        if se == 'Start':
            tripLocation = 11  
            tripTime = 3
            tripHomeID = 1
            tripDist = 5
            TripLatitude = 7
            TripLongitude = 8 
            TripClusterID = 13
        else:
            tripLocation = 12
            tripTime = 4
            tripHomeID = 2
            tripDist = 6
            TripLatitude = 9
            TripLongitude = 10
            TripClusterID = 14
        for row in vLocationInfo.itertuples():
            vrow = vLocationInfo.iloc[row.Index]

            ### Set Trip(Start/End)LocationCategory and Trip(Start/End)HomeID
            # Does row match match a HomeID and is it in between any Home(Start/Local)LocalTime ranges
            homeID, inRange = isInHomeCluster(se, row[TripClusterID+1], row[tripTime+1], homeInfo)
            ## did not match a home
            if homeID == -1:
                category = 'unknown'
                # is park in any home cluster
                if isInRangeSet(se, homeInfo[['HomeStartLocalTime', 'HomeEndLocalTime']], row[tripTime+1]):
                    category = 'away'
                vrow[tripLocation] = category
                vrow[tripHomeID] = None

            ## matched a home and is in home cluster period
            if homeID != -1 and inRange:
                vrow[tripLocation] = 'home'
                vrow[tripHomeID] = homeID

            ## matched a home, but is not in home cluster period
            if homeID != -1 and inRange == False:
                category = 'unknown'
                # is park in any home cluster
                if isInRangeSet(se, homeInfo[['HomeStartLocalTime', 'HomeEndLocalTime']], row[tripTime+1]):
                    category = 'away'
                vrow[tripLocation] = category
                vrow[tripHomeID] = None

            homeLoc = []
            for homerow in homeInfo[homeInfo['Primary'] == True].itertuples():
                if isInTupleRange(se, homerow.HomeStartLocalTime, homerow.HomeEndLocalTime, row[tripTime+1]):
                    homeLoc.extend([homerow[0]])

            # no homes
            if len(homeLoc) == 0:
                vrow[tripDist] = None
            # one distinct home
            if len(homeLoc) == 1:
                if row[tripLocation] == 'home':
                    vrow[tripDist] = 0
                else:
                    hm = homeInfo.iloc[homeLoc[0]]
                    vrow[tripDist] = math.ceil(haversine((row[TripLatitude+1], row[TripLongitude+1]), (hm['CentroidLatitude'], hm['CentroidLongitude']), unit=Unit.MILES))
            # in range with multiple homes
            if len(homeLoc) > 1:
                hm = homeInfo[homeInfo.index.isin(homeLoc) & (homeInfo['Primary'])]
                vrow[tripDist] = math.ceil(haversine((row[TripLatitude+1], row[TripLongitude+1]), (hm['CentroidLatitude'], hm['CentroidLongitude']), unit=Unit.MILES))

            vLocationInfo.iloc[row.Index] = vrow
    return True, vLocationInfo, homeInfo

def flagTrips(v, vData):
    # use offset as end/start of day, e.g. 3:30 AM
    vData['TripStartDateOffset'] = (vData['TripStartLocalTime'] - timedelta(hours=cfg.dayEndHours, minutes=cfg.dayEndMinutes)).dt.date
    vData['TripEndDateOffset']= (vData['TripEndLocalTime'] - timedelta(hours=cfg.dayEndHours, minutes=cfg.dayEndMinutes)).dt.date

    lastIdx = len(vData) - 1
    curParkEndDate = vData['TripStartDateOffset'][0:1].item()
    vData['TripFlag'][0:1] = 'F'
    tripsCnt = 0
    # find first and last trips in the day
    for i in range(1, lastIdx):
        tripsCnt += 1
        # compare current (i) record to endDate
        if vData['TripEndDateOffset'][i:i+1].item() != curParkEndDate:   
            vData['TripFlag'][i-1:i] = 'FL' if vData['TripFlag'][i-1:i].item() == 'F' else 'L' 
            vData['TripFlag'][i:i+1] = 'F'
            curParkEndDate = vData['TripEndDateOffset'][i:i+1].item()
            tripsCnt = 0
    vData['TripFlag'][-1:] = 'FL' if vData['TripFlag'][lastIdx-1:lastIdx].item() == 'L' else 'L'

    return vData

def InBoundingBox(vd, colLat, colLon):
    """Check a value (latitude or longitude) to see if it is within the given range"""
    if math.isnan(vd[colLat]) or math.isnan(vd[colLon]):
        vd['Include'] = False
        return vd

    x = vd[colLat]
    y = vd[colLon]
    isFound = False
    for k in cfg.bboxes.keys():
        x1 = cfg.bboxes[k][0][0] # upper-
        y1 = cfg.bboxes[k][0][1] #   left coordinates
        x2 = cfg.bboxes[k][1][0] # lower- 
        y2 = cfg.bboxes[k][1][1] #   right coordinates
        if x > x2 and x < x1 and y > y1 and y < y2: # note: x-axis decreases from bottom to top
            isFound = True
            break

    # don't change any previously "falsed" flags
    if not isFound:
        vd['Include'] = False
    
    return vd

# check that dates and times are sane
def CheckDateTime(vd, colname):
    try:
        if pd.isnull(vd[colname]):
            vd['Include'] = False
            return vd
        curdt = datetime.today()
        if vd[colname].year < 2011 or vd[colname] > curdt:
            vd['Include'] = False
        return vd
    except ValueError:
        vd['Include'] = False
        return vd

# check that the Odometer mileage is not less than the calculated mileage from the GPS coordinates
def CompareOdometerToGPS(vd, tripStart, tripEnd, stLat, stLng, enLat, enLng, threshold):
    odoDist = vd[tripEnd] - vd[tripStart]
    GPSDist = haversine((vd[stLat], vd[stLng]), (vd[enLat], vd[enLng]), unit=Unit.MILES)
    if (odoDist - GPSDist) < threshold:
        vd.Include = False
    return vd

# check various types of data with the vehicle data data frame and return data frame with the Include flag set
def DoValidityChecking(v, vData):
    incl = vData

    incl = incl.apply(lambda x: CheckDateTime(x, 'TripStartLocalTime'), axis=1)
    startErrs = incl['Include'][incl['Include'] == False].count()
    if startErrs > 0: cfg.errorWriter.writerow([datetime.now(), v, -1, 'TripStartLocalTimes ({})'.format(startErrs)])

    incl = incl.apply(lambda x: CheckDateTime(x, 'TripEndLocalTime'),   axis=1)
    endErrs = incl['Include'][incl['Include'] == False].count() - startErrs
    if endErrs > 0: cfg.errorWriter.writerow([datetime.now(), v, -1, 'TripEndLocalTime is bad ({})'.format(endErrs)])

    incl = incl.apply(lambda x: InBoundingBox(x, 'TripStartLatitude', 'TripStartLongitude'), axis=1)
    startPtErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs
    if startPtErrs > 0: cfg.errorWriter.writerow([datetime.now(), v, -1, 'TripStartLatitude is bad ({})'.format(startPtErrs)])

    incl = incl.apply(lambda x: InBoundingBox(x, 'TripEndLatitude', 'TripEndLongitude'), axis=1)
    endPtErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs - startPtErrs
    if endPtErrs > 0: cfg.errorWriter.writerow([datetime.now(), v, -1, 'TripEndLongitude is bad ({})'.format(endPtErrs)])

    incl = incl.apply(lambda x: CompareOdometerToGPS(x, 'TripStartOdometer_Miles', 'TripEndOdometer_Miles', 
                                                         'TripStartLatitude', 'TripStartLongitude', 
                                                         'TripEndLatitude', 'TripEndLongitude', cfg.gpsOdoThreshold_mi), axis=1)
    odoErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs - startPtErrs - endPtErrs
    if odoErrs > 0: cfg.errorWriter.writerow([datetime.now(), v, -1, 'Trip ODO < straight line distance ({})'.format(odoErrs)])
    
    return incl

# find clusters of data using latitude and longitude of trip start and end
def clusterData(v,vData):
    kms_per_radian = 6371.0088
    epsilon = (cfg.dbscan_eps_ft / 3281) / kms_per_radian
    minSamples = cfg.dbscan_min_spls

    startPts = vData[['TripStartLatitude', 'TripStartLongitude']].to_numpy()
    endPts = vData[['TripEndLatitude', 'TripEndLongitude']].to_numpy()
    coords = np.append(startPts, endPts, axis=0)
    coordsSet = np.unique(coords, axis=0)

    db = DBSCAN(eps=epsilon, min_samples=minSamples, algorithm='ball_tree', metric='haversine').fit(np.radians(coordsSet))
    clusterLbls = db.labels_ # db.labels seems to be an array of cluster IDs mapping to the coords array

    coordsClusters = pd.DataFrame(coordsSet)
    coordsClusters['clusterID'] = clusterLbls
    coordsClusters.columns = ['latitude', 'longitude', 'clusterID']

    # map vData locations to their cluster ID
    vData['TripStartClusterID'] = None
    vData['TripEndClusterID'] = None
    for vIdx, row in vData.iterrows():
        cIdx = coordsClusters.index[(coordsClusters['latitude'] == row['TripStartLatitude']) & (coordsClusters['longitude'] == row['TripStartLongitude'])]
        vData['TripStartClusterID'][vIdx] = coordsClusters['clusterID'][cIdx].values[0]
        cIdx = coordsClusters.index[(coordsClusters['latitude'] == row['TripEndLatitude']) & (coordsClusters['longitude'] == row['TripEndLongitude'])]
        vData['TripEndClusterID'][vIdx] = coordsClusters['clusterID'][cIdx].values[0]

    return vData

if __name__ == '__main__':

    main()
