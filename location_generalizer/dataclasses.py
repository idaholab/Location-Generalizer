from dataclasses import dataclass

@dataclass
class StartSEColumnMappings:
    '''Column names for a start SE'''
    mapType: str = 'Start'
    tripLocation: str = 'TripStartLocationCategory'  
    tripTime: str = 'TripStartLocalTime'
    tripHomeID : str = 'TripStartHomeID'
    tripDistance : str = 'TripStartDistanceFromHome_Miles'
    tripLatitude : str = 'TripStartLatitude'
    tripLongitude : str = 'TripStartLongitude'
    tripClusterID : str = 'TripStartClusterID'

@dataclass
class EndSEColumnMappings:
    '''Column names for an end SE'''
    mapType: str = 'End'
    tripLocation: str = 'TripEndLocationCategory'  
    tripTime: str = 'TripEndLocalTime'
    tripHomeID : str = 'TripEndHomeID'
    tripDistance : str = 'TripEndDistanceFromHome_Miles'
    tripLatitude : str = 'TripEndLatitude'
    tripLongitude : str = 'TripEndLongitude'
    tripClusterID : str = 'TripEndClusterID'
