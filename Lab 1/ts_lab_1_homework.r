###################################################################
# preliminary analysis of ariline passengers time series 
###################################################################

# add a comment after every # in some cases there are indications 

library(TSA)
data()

data(AirPassengers)
y = AirPassengers;

# describe the data 
y


# preliminary analysis of the graph
ts.plot(y)

# add an extensive comment here on the choice of the log-transform
# you can do a range-mean plot if you want 

ly = log(y)

# 
ts.plot(ly)

#
acf(ly, lag.max = 48)

#
dly = diff(ly)

#
ts.plot(dly)

#
acf(dly, lag.max = 48)


#  
ddly = diff(dly,12)
ts.plot(ddly)

#  
dsly = diff(ly,12)
ts.plot(dsly)
acf(dsly, lag.max = 11)

#
ddsly = diff(dsly)
ts.plot(ddsly)




