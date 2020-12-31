########################################################################
#  Copyright 2020 Battelle Energy Alliance, LLC   ALL RIGHTS RESERVED  #
#  Mobility Systems & Analytics Group, Idaho National Laboratory       #
########################################################################

import pyodbc
import pandas as pd
import pickle
import datetime
import time
import math
import yaml
#import geopandas
#import shapely
from pathlib import Path
import csv
import numpy as np
from sklearn.cluster import DBSCAN
from shapely import geometry
from shapely.geometry import MultiPoint
from haversine import haversine, Unit
import pynput

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
    distancePlaces = config['distancePlaces']

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

    verbose = 0
    stopProcessing = False

    errFilePath = Path(errorLogFileName)
    if not errFilePath.exists():
        # ErroLog output file
        hdr = pd.DataFrame(hdrErrorLogCSV)
        hdr.to_csv(errorLogFileName, index=False, header=False, mode='w')
    # use one line buffering - every line written is flushed to disk
    errorFile = open(errorLogFileName, mode='a', buffering=1, newline='')
    errorWriter = csv.writer(errorFile)

def on_press(key):
    if hasattr(key, 'char'):
        if key.char == 'v':
            cfg.verbose = (cfg.verbose + 1) % 3 # verbosity levels: 0, 1, 2
            print('Verbosity: {}'.format(cfg.verbose))
        if key.char == 'q':
            cfg.stopProcessing = not cfg.stopProcessing
            if cfg.stopProcessing:
                print("Processing will stop after current vehicle.")
            else:
                print("Stop canceled, processing will continue.")

def main():
    listener = pynput.keyboard.Listener(on_press=on_press)
    listener.start() 

    # for vehicle proceesing rate
    vst = datetime.datetime.now()

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

    # get the census division polygons to use with finding home
    # divisions = geopandas.read_file(cfg.censusDivisionsFileName)
    # divisions = geopandas.GeoDataFrame.from_file(cfg.shpCensusDivisionsFileName)
    # divisions.to_pickle(cfg.pklCensusDivisionsFileName)

    ## geopandas can read the shapefile directly, but we pickled it into one file
    ## a single pickle file simplifies distribution whereas,
    ## loading a shapefile requires several adjacent accompanying files
    divisions = pd.read_pickle(cfg.pklCensusDivisionsFileName)

    # get Public EVSE stations
    EVSEs = pd.read_csv(cfg.evseLookupFileName)

    # create empty LocationInfo data frame
    # GPS coordinates are added here for convenience, but will not be carried into LocationInfo output file
    locationInfo = pd.DataFrame(columns = cfg.colLocationInfo)
    
    # create empty HomeInfo data frame
    homeInfo = pd.DataFrame(columns = cfg.colHomeInfo)

    # autocommit here is a workaround to prevent pyodbc from erroring when connecting to CSV files
    pyodbc.autocommit = True
    cnxn = pyodbc.connect(cfg.odbcConnectionString, autocommit=True)


    lastVehicle = 0
    hbFilePath = Path(cfg.heartbeatFileName)
    if hbFilePath.exists():
        with open(hbFilePath, 'r') as hb:
            lastVehicle = hb.readline()
        cfg.errorWriter.writerow([datetime.datetime.now(), lastVehicle, -1,'Restarting after vehicle {}'.format(lastVehicle)])
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
        if cfg.verbose > 0: print('Fetching data')
        chunkData = pd.read_sql(qry, cnxn, parse_dates=['TripStartLocalTime', 'TripEndLocalTime'])
        # sqlQry = pd.read_sql_query(qry, cnxn)
        # chunkData = pd.DataFrame(sqlQry)
        # create new column for flag to exclude bad records
        chunkData['Include'] = True
        i += 1
        print("chunk: {}, vehicle from {} through {}".format(i, chunk[0], chunk[-1]))

        # iterate through one vehicle at a time
        for v in chunk:
            if cfg.stopProcessing: exit()
            if cfg.verbose > 0: print('Vehicle: {}'.format(v))
            vcnt += 1
            # grab all records in vehicle v
            vData = chunkData[chunkData['VehicleID'] == v]
            # create new column to check for Odometer gaps, i.e missing trips
            vData['resid_Miles'] = vData['TripStartOdometer_Miles'].shift(periods=-1) - vData['TripEndOdometer_Miles']

            ### Check validity of data, marking invalid records (Include = True/False)
            if cfg.verbose > 1: print(' Check for valid values')
            vData = DoValidityChecking(v, vData)
            vData.resid_Miles = vData.resid_Miles.astype(object).where(vData.resid_Miles.notnull(), None) # set NaN to None (becomes Null for DB)
            # toss out rows that failed vailidity check
            vData = vData[vData.Include == True]
            numTrips = len(vData)
            if numTrips < cfg.minTrips:
                if cfg.verbose > 1: print(' Not enough trips, vehicle skipped.')
                cfg.errorWriter.writerow([datetime.datetime.now(), v, -1,'Not enough trips, vehicle skipped. ({} need >= {}'.format(numTrips, cfg.minTrips)])
            else:
                # create new column for identify first/last trip of day
                vData['TripFlag'] = None
                ### Identify first and last of trip of day
                if cfg.verbose > 1: print(' Defining first/last trip of day')
                vData = flagTrips(v, vData)
                ### Find clusters of vehicle locations
                if cfg.verbose > 1: print(' Clustering')
                vData, clusterData(v, vData)

                # drop rows - remove previous vehicle info
                homeInfo.drop(homeInfo.index, inplace=True)
                locationInfo.drop(locationInfo.index, inplace=True)

                # add row to LocationInfo data frame
                liList = [vData[['VehicleID', 'TripStartLocalTime', 'TripEndLocalTime', 'TripStartLatitude', 'TripStartLongitude', 'TripEndLatitude','TripEndLongitude', 'TripStartClusterID', 'TripEndClusterID']]]
                locationInfo = locationInfo.append(liList, ignore_index=True)

                ########################
                #### FIND HOME MODULE 1: must have at least 30 valid last trip of day and 95% of those in one cluster
                #### (does not actually return more than one home, but must returns an array to conform with other methods that may return more that one)
                if cfg.verbose > 1: print(' Identifying home location')
                topClusterIDs, homeCPs = findHome_Module1(v, vData)
                ########################

                # process location and home info for home locations
                if cfg.verbose > 1: print(' Calculating output metrics')
                locationInfo, homeInfo = processHome(v, divisions, vData, locationInfo, homeInfo, topClusterIDs, homeCPs, EVSEs)

                if not homeInfo.empty:
                    # write to output files
                    if cfg.verbose > 1: print(' Writing to output files')
                    locationInfo.to_csv(cfg.locationInfoFileName, index=False, header=False, mode='a')
                    homeInfo.to_csv(cfg.homeInfoFileName, index=False, header=False, mode='a')
            # # use one line buffering - every line written is flushed to disk
            with open(cfg.heartbeatFileName, mode='w', buffering=1, newline='') as hb:
                hb.write(str(v))
            if vcnt % 5 == 0:
                ven = datetime.datetime.now()
                secs = (ven - vst).seconds
                if cfg.verbose > 0: 
                    if secs > 0:
                        rate = 3600 * (vcnt/secs)
                        print('Processing rate (cumulative average): {:4.0f} / hour   '.format(rate))

def findHome_Module1(v, vData):
    CP = geometry.Point(0,0)
    topClusterID = -1
    # get the parks after known last trip of day
    lastParks = vData[((vData['TripFlag'] == 'L') | (vData['TripFlag'] == 'FL')) & (vData['resid_Miles'] < 1) & (vData['resid_Miles'] > -1)]
    if not lastParks.empty:
        numValidLastTrips = lastParks['TripFlag'].count()
        # get top cluster with last parks
        #z = lastParks.groupby('TripEndClusterID')['TripEndClusterID'].count()
        z = lastParks.groupby('TripEndClusterID')['TripEndClusterID'].count().sort_values(ascending=False)
        topClusterCnt = z.iloc[0]
        topClusterID = z.index[0]

        if numValidLastTrips < cfg.minLastTrips:
            cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'Not enough last trips ({})'.format(numValidLastTrips)])
        else:
            #if topClusterCnt <= (numValidLastTrips *  cfg.minPctParks):
                # cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'No home by percent trips. {} in clusterID {}, needed {:1.0f}% of {} ({:1.0f}).'.format(topClusterCnt, topClusterID, cfg.minPctParks*100,  numValidLastTrips,  (numValidLastTrips *  cfg.minPctParks))])
            #else:
            if topClusterCnt > (numValidLastTrips *  cfg.minPctParks):
                # get centroid of top cluster (use only the home filtered points)
                ## get the home filtered points belonging to the topClusterID
                dfPts = lastParks[['TripEndLatitude', 'TripEndLongitude']][lastParks['TripEndClusterID'] == topClusterID]
                ## convert to MultiPoint object
                mpPts = MultiPoint(dfPts.to_numpy())
                ## get the center point of the topCLusterID
                CP = mpPts.centroid
                # need lat/long positions switched for "within" check
                CP = geometry.Point(CP.y, CP.x)
    # return as arrays to conform with other findHome modules
    homeCPs = [CP, geometry.Point(0,0)]
    topClusterIDs = [topClusterID, -1]
    return topClusterIDs, homeCPs

# moving window module
def findHome_Module2(v, vData):
    print(vData)
    for i in range(0, len(vData), 30): print(i)

def getEVSEDistance(row, homeLat, homeLong):
    dist = haversine((row.Latitude, row.Longitude), (homeLat, homeLong), unit=Unit.MILES)
    return dist

def getStartLocationDistance(row, homeLat, homeLong, homeStart, homeEnd):
    if (homeStart <= row['TripStartLocalTime'] <= homeEnd):
        startDist = haversine((row['TripStartLatitude'], row['TripStartLongitude']), (homeLat, homeLong), unit=Unit.MILES)
    else:
        startDist = row['TripStartDistanceFromHome_Miles']
    return startDist

def getEndLocationDistance(row, homeLat, homeLong, homeStart, homeEnd):
    endDist = haversine((row['TripEndLatitude'], row['TripEndLongitude']), (homeLat, homeLong), unit=Unit.MILES)
    if (homeStart <= row['TripEndLocalTime'] <= homeEnd):
        endDist = haversine((row['TripEndLatitude'], row['TripEndLongitude']), (homeLat, homeLong), unit=Unit.MILES)
    else:
        endDist = row['TripEndDistanceFromHome_Miles']
    return endDist

def processHome(v, divisions, vData, vLocationInfo, homeInfo, topClusterIDs, homeCPs, EVSEs):
    # loop through home(s) from current vehicle
    for cIdx, homeCPgeom in enumerate(homeCPs):
        if homeCPgeom.x == 0.0: continue
        #homeCPgeom = geometry.Point(CP[1], CP[0])
        for i, division in divisions.iterrows():
            if division.geometry.contains(homeCPgeom):
                st = EVSEs[(EVSEs['Latitude'] > (homeCPgeom.y - cfg.evseLatRange)) & 
                            (EVSEs['Latitude'] < (homeCPgeom.y + cfg.evseLatRange)) & 
                            (EVSEs['Longitude'] > (homeCPgeom.x - cfg.evseLonRange)) & 
                            (EVSEs['Longitude'] < (homeCPgeom.x + cfg.evseLonRange))]
                if not st.empty:
                    st['hMiles'] = st.apply(getEVSEDistance, args=(homeCPgeom.y, homeCPgeom.x), axis=1)
                    st = st[st['hMiles'] <= cfg.evseDistRange_Miles]
                l2Cnt = 0
                dcCnt = 0
                if not st.empty:
                    l2Cnt = st['L2'].sum()
                    dcCnt = st['DCFC'].sum()

                # add info to homeInfo
                #homeStart = vData['TripStartLocalTime'][vData['ClusterID'] == topClusterIDs[cIdx]].min()
                #homeEnd = vData['TripEndLocalTime'][vData['ClusterID'] == topClusterIDs[cIdx]].max()

                homeStart = vData['TripStartLocalTime'].min()
                homeEnd = vData['TripEndLocalTime'].max()
                newRow = {'VehicleID':v, 
                          'HomeStartLocalTime':homeStart, 'HomeEndLocalTime':homeEnd, 
                          'HomeRegion':division['NAME'], 'PublicChargingDensityL2':l2Cnt, 'PublicChargingDensityDCFC':dcCnt}
                homeInfo = homeInfo.append(newRow, ignore_index=True)

                # compute distance from home to trip start/end
                vLocationInfo['TripStartDistanceFromHome_Miles'] = vLocationInfo.apply(getStartLocationDistance, args=(homeCPgeom.y, homeCPgeom.x, homeStart, homeEnd), axis = 1)
                vLocationInfo['TripEndDistanceFromHome_Miles'] = vLocationInfo.apply(getEndLocationDistance, args=(homeCPgeom.y, homeCPgeom.x, homeStart, homeEnd), axis = 1)
                vLocationInfo = vLocationInfo.round({'TripStartDistanceFromHome_Miles': cfg.distancePlaces, 'TripEndDistanceFromHome_Miles': cfg.distancePlaces})

                # categorize locations
                vLocationInfo.loc[vLocationInfo['TripStartClusterID'] == topClusterIDs[cIdx], 'TripStartLocationCategory'] = 'home'
                vLocationInfo.loc[vLocationInfo['TripEndClusterID'] == topClusterIDs[cIdx], 'TripEndLocationCategory'] = 'home'
                vLocationInfo.loc[vLocationInfo['TripStartClusterID'] != topClusterIDs[cIdx], 'TripStartLocationCategory'] = 'away'
                vLocationInfo.loc[vLocationInfo['TripEndClusterID'] != topClusterIDs[cIdx], 'TripEndLocationCategory'] = 'away'
               
                break

        # remove GPS location info
        vLocationInfo.drop(['TripStartLatitude', 'TripStartLongitude', 'TripEndLatitude', 'TripEndLongitude'], axis=1, inplace=True)

    return vLocationInfo, homeInfo

def flagTrips(v, vData):
    # use offset as end/start of day, e.g. 3:30 AM
    vData['TripStartDateOffset'] = (vData['TripStartLocalTime'] - datetime.timedelta(hours=cfg.dayEndHours, minutes=cfg.dayEndMinutes)).dt.date
    vData['TripEndDateOffset']= (vData['TripEndLocalTime'] - datetime.timedelta(hours=cfg.dayEndHours, minutes=cfg.dayEndMinutes)).dt.date

    lastIdx = len(vData) - 1
    curParkEndDate = vData['TripStartDateOffset'][0:1].item()
    vData['TripFlag'][0:1] = 'F'
    tripsCnt = 0
    # find first and last trips in the day
    for i in range(1, lastIdx):
        tripsCnt += 1
        # compare current (i) record to endDate
        if vData['TripEndDateOffset'][i:i+1].item() != curParkEndDate:   
            vData['TripFlag'][i-1:i] = 'FL' if vData['TripFlag'][i-1:i].item() == 'F' else 'L'    # 
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

def CheckDateTime(vd, colname):
    try:
        if pd.isnull(vd[colname]):
            vd['Include'] = False
            return vd
        curdt = datetime.datetime.today()
        if vd[colname].year < 2011 or vd[colname] > curdt:
            vd['Include'] = False
        return vd
    except ValueError:
        vd['Include'] = False
        return vd

def CompareOdometerToGPS(vd, tripStart, tripEnd, stLat, stLng, enLat, enLng, threshold):
    """Check that the Odometer mileage is not less than the calculated mileage from the GPS coordinates"""
    odoDist = vd[tripEnd] - vd[tripStart]
    GPSDist = haversine((vd[stLat], vd[stLng]), (vd[enLat], vd[enLng]), unit=Unit.MILES)
    if (odoDist - GPSDist) < threshold:
        vd.Include = False
    return vd

def DoValidityChecking(v, vData):
    """Check various types of data with the vehicle data data frame and return data frame with the Include flag set"""
    incl = vData

    incl = incl.apply(lambda x: CheckDateTime(x, 'TripStartLocalTime'), axis=1)
    startErrs = incl['Include'][incl['Include'] == False].count()
    if startErrs > 0: cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'TripStartLocalTimes ({})'.format(startErrs)])

    incl = incl.apply(lambda x: CheckDateTime(x, 'TripEndLocalTime'),   axis=1)
    endErrs = incl['Include'][incl['Include'] == False].count() - startErrs
    if endErrs > 0: cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'TripEndLocalTime is bad ({})'.format(endErrs)])

    incl = incl.apply(lambda x: InBoundingBox(x, 'TripStartLatitude', 'TripStartLongitude'), axis=1)
    startPtErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs
    if startPtErrs > 0: cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'TripStartLatitude is bad ({})'.format(startPtErrs)])

    incl = incl.apply(lambda x: InBoundingBox(x, 'TripEndLatitude', 'TripEndLongitude'), axis=1)
    endPtErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs - startPtErrs
    if endPtErrs > 0: cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'TripEndLongitude is bad ({})'.format(endPtErrs)])

    incl = incl.apply(lambda x: CompareOdometerToGPS(x, 'TripStartOdometer_Miles', 'TripEndOdometer_Miles', 
                                                         'TripStartLatitude', 'TripStartLongitude', 
                                                         'TripEndLatitude', 'TripEndLongitude', cfg.gpsOdoThreshold_mi), axis=1)
    odoErrs = incl['Include'][incl['Include'] == False].count() - startErrs - endErrs - startPtErrs - endPtErrs
    if odoErrs > 0: cfg.errorWriter.writerow([datetime.datetime.now(), v, -1, 'Trip ODO < straight line distance ({})'.format(odoErrs)])
    
    return incl

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
