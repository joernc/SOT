using CSV, DataFrames

# load ISC catalog, select required columns
iscdf = DataFrame(CSV.File("catalogs/ISC_nias.csv", header=26, delim=",",
                           select=[1, 3, 4, 5, 6], normalizenames=true))

# generate new DataFrame, convert data to desired format
newdf = DataFrame()
newdf.eventid = iscdf.EVENTID
newdf.time = iscdf.DATE + iscdf.TIME
newdf.latitude = iscdf.LAT
newdf.longitude = iscdf.LON

# write to file
CSV.write("catalogs/nias.csv", newdf)
