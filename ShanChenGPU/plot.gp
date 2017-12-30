set term gif
set output "$0.gif"
set zrange [0.0:2.3]
splot [0:255] [0:255] '$0' matrix w l
#!kuickshow test.jpg
reset