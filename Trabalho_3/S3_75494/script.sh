g++ S3_onda.cpp -o S3_onda
./S3_onda

#Output splot
#gnuplot Output_splot.plt

#Output gif
gnuplot Output_gif.gp
xdg-open Output_gif.gif
#se quisermos aumentar a velocidade do gif
#convert -delay 1x30 foobar.gif mod.gif 
