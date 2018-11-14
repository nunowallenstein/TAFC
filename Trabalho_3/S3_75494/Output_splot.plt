
#set zrange[0:1.]
#set title "Gráfico da onda na equação de Maxwell,\n lambda =0 (Sobreposição de Modos N_t=100000) \n Condições fronteira reflectoras"
set xlabel "Position"
set ylabel "Time"
#set zlabel "Psi"
set ticslevel 0
set pm3d
set palette defined (-3 "white",0 "green", 1 "red")
unset border;
unset ztics;
#gráfico a plotar

#right(lambda =-1)
#periodicas 	 
#set view 11,68,1,1
#reflectoras
set view 60,80,1,1


#central periodicas
#set view 11,311,1,1
#central reflectoras
#set view 75,20,1

#left periodicas
#set view 30,280,1,1
#reflectoras
#set view 30,320,1,1

splot 'Output_psi.dat' matrix nonuniform with dots notitle

	
pause -1 "Hit return to exit"
