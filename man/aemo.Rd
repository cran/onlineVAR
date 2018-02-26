\name{aemo}
\alias{aemo}

\title{Wind power data from 22 wind farms in south eastern Australia}

\description{
Wind power data from 22 wind farms in south eastern Australia provided by
the Australian Electricity Market Operator (AEMO) and pre-processed by Jethro 
Browell and Stefanos Delikaraoglou. The data set contains data from 2013
in 5 minute resolution and is a subset of the data set provided 
at http://dx.doi.org/10.15129/9e1d9b96-baa7-4f05-93bd-99c5ae50b141.
The columns are the wind farms sorted from west to east.
Meta data on the capacity and location of the wind farms can be accessed 
with \code{attr(aemo, "meta")}.
}

\usage{data("aemo")}

\format{ 
  A data frame with 22 columns for each wind farm and 105120 rows. 
}

\source{
http://dx.doi.org/10.15129/9e1d9b96-baa7-4f05-93bd-99c5ae50b141
https://benjaminjweise.carto.com/tables/aemo_wind_plants/public
}

\references{
  Dowell J, Pinson P (2016), Very-Short-Term Probabilistic Wind Power Forecasts  
  by Sparse Vector Autoregression. \emph{IEEE Transactions on Smart Grid}, 7(2),
  763-770, DOI: 10.1109/TSG.2015.2424078
}

\keyword{datasets}
