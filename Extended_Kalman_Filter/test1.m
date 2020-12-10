directory="D:\Documents\Thesis_Research\MIT_Sailing\processed_data\quokka";
start_time=1536932173.36;
s=geoshape(NAVLL.NAV_LAT(5553:37948),NAVLL.NAV_LONG(5553:37948));
filename_kml=directory + string(start_time) + "_NAV1.kml";
kmlwrite(filename_kml,s);