library(rdetools)

?sinc
m = 1/2
h = 2*3.47088972434702
zero = 0.0913252
phi01 = function(x) sinc((x-m)*h) + zero #function(x) sin(x)/x - 7.72525183693771
plot(phi01, from=0, to =1)


norm_phi01 = sqrt(integrate(function(x){phi01(x)^2}, lower = 0, upper = 1)$value)

phi1 = function(x) phi01(x)/norm_phi01


phi2 = function(x) sqrt(2)*sin(8*pi*x)
plot(phi2, from=0,to=1)
