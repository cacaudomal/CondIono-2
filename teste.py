# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 16:44:08 2023

@author: tedea
"""

import matplotlib.pyplot as plt
import pandas as pd
import pyigrf_clara_0_4 as igrf
import condutividade_0_7 as cd

name_saida = "teste_cond"
hour = ["00LT","06LT","12LT","18LT"]
month = ['JAN',"MAR","SEP",'DEZ']

IRIpath2000 = ["IRI2016_lat-25_lon125_z0_745_s_5data01JAN2000_h" + i +".txt" for i in hour]
MSISEpath2000 = ["NRLMSISE2_lat-25_lon125_z0_745_s_5data01JAN2000_h" + i + ".txt" for i in hour]

IRIpath2008 = ["IRI2016_lat-25_lon125_z0_745_s_5data01JAN2008_h" + i +".txt" for i in hour]
MSISEpath2008 = ["NRLMSISE2_lat-25_lon125_z0_745_s_5data01JAN2008_h" + i + ".txt" for i in hour]

IRIpath2008mes = ["IRI2016_lat-25_lon125_z0_745_s_5data01" + i + "2008_h00LT.txt" for i in month]
MSISEpath2008mes = ["NRLMSISE2_lat-25_lon125_z0_745_s_5data01" + i + "2008_h00LT.txt" for i in month]

igrf_2000 = igrf.IGRF(-25,55,0,2000,name_saida)
igrf_2000.calc_perfilh(-25,55,lim_h = 750)

igrf_2008 = igrf.IGRF(-25,-55,0,2000,name_saida)
igrf_2008.calc_perfilh(-25,-55,lim_h = 750)

cond2000 = [cd.condiono_adachi(igrf_2000.DfPerfilh["B(T)"],IRIpath2000[i],MSISEpath2000[i]) for i in range(len(hour))]
cond2008 = [cd.condiono_adachi(igrf_2008.DfPerfilh["B(T)"],IRIpath2008[i],MSISEpath2008[i]) for i in range(len(hour))]
cond2008mes = [cd.condiono_adachi(igrf_2008.DfPerfilh["B(T)"],IRIpath2008mes[i],MSISEpath2008mes[i]) for i in range(len(month))]


condHall2000 = [ cond2000[i].calc_Hall(cond2000[i].Freq["fen"],cond2000[i].Freq["fin1"],cond2000[i].Freq["fin2"],\
                    cond2000[i].girofreq["wi1"],cond2000[i].girofreq["wi2"],cond2000[i].girofreq["we"],
                    cond2000[i].p1, cond2000[i].p2, cond2000[i].iri.dado["Ne(m-3)"], cond2000[i].B)\
                for i in range(len(hour))]
condPedersen2000 = [ cond2000[i].calc_Pedersen(cond2000[i].Freq["fen"],cond2000[i].Freq["fin1"],cond2000[i].Freq["fin2"],\
                  cond2000[i].girofreq["wi1"],cond2000[i].girofreq["wi2"],cond2000[i].girofreq["we"],\
                  cond2000[i].p1, cond2000[i].p2, cond2000[i].iri.dado["Ne(m-3)"], cond2000[i].B)
                for i in range(len(hour))]


condHall2008 = [ cond2008[i].calc_Hall(cond2008[i].Freq["fen"],cond2008[i].Freq["fin1"],cond2008[i].Freq["fin2"],\
                    cond2008[i].girofreq["wi1"],cond2008[i].girofreq["wi2"],cond2008[i].girofreq["we"],
                    cond2008[i].p1, cond2008[i].p2, cond2008[i].iri.dado["Ne(m-3)"], cond2008[i].B)\
                for i in range(len(hour))]
condPedersen2008 = [ cond2008[i].calc_Pedersen(cond2008[i].Freq["fen"],cond2008[i].Freq["fin1"],cond2008[i].Freq["fin2"],\
                  cond2008[i].girofreq["wi1"],cond2008[i].girofreq["wi2"],cond2008[i].girofreq["we"],\
                  cond2008[i].p1, cond2008[i].p2, cond2008[i].iri.dado["Ne(m-3)"], cond2008[i].B)
                for i in range(len(hour))]


condHall2008mes = [ cond2008mes[i].calc_Hall(cond2008mes[i].Freq["fen"],cond2008mes[i].Freq["fin1"],cond2008mes[i].Freq["fin2"],\
                    cond2008mes[i].girofreq["wi1"],cond2008mes[i].girofreq["wi2"],cond2008mes[i].girofreq["we"],
                    cond2008mes[i].p1, cond2008mes[i].p2, cond2008mes[i].iri.dado["Ne(m-3)"], cond2008mes[i].B)\
                for i in range(len(hour))]    
condPedersen2008mes = [ cond2008mes[i].calc_Pedersen(cond2008mes[i].Freq["fen"],cond2008mes[i].Freq["fin1"],cond2008mes[i].Freq["fin2"],\
                  cond2008mes[i].girofreq["wi1"],cond2008mes[i].girofreq["wi2"],cond2008mes[i].girofreq["we"],\
                  cond2008mes[i].p1, cond2008mes[i].p2, cond2008mes[i].iri.dado["Ne(m-3)"], cond2008mes[i].B)
                for i in range(len(hour))]

    
# fig4,(ax1,ax2) = plt.subplots(1,2,figsize=(10, 5))
# for i in range(len(hour)):
#     ax1.semilogx(condHall2008mes[i]*-1,cond2008mes[i].iri.dado["H(km)"],label = month[i])
# ax1.legend() 
# ax1.set_title("Hall")

# for i in range(len(hour)):
#     ax2.semilogx(condPedersen2008mes[i],cond2008mes[i].iri.dado["H(km)"],label = month[i])
# ax2.legend() 
# ax2.set_title("Pedersen")
# fig4.suptitle("Conductivities month \nyear: 2008 lat: -25 lon: 125")
    
    
# fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10, 5))
# for i in range(len(hour)):
#     ax1.semilogx(condHall2000[i]*-1,cond2000[i].iri.dado["H(km)"],label = hour[i])
# ax1.legend() 
# ax1.set_title("Hall")

# for i in range(len(hour)):
#     ax2.semilogx(condPedersen2000[i],cond2000[i].iri.dado["H(km)"],label = hour[i])
# ax2.legend() 
# ax2.set_title("Pedersen")
# fig.suptitle("Conductivities \nyear: 2000 lat: -25 lon: 125")

#==== 
# fig2,(ax1,ax2) = plt.subplots(1,2,figsize=(10, 5))
# for i in range(len(hour)):
#     ax1.semilogx(condHall2008[i]*-1,cond2008[i].iri.dado["H(km)"],label = hour[i])
# ax1.legend() 
# ax1.set_title("Hall")

# for i in range(len(hour)):
#     ax2.semilogx(condPedersen2008[i],cond2008[i].iri.dado["H(km)"],label = hour[i])
# ax2.legend() 
# ax2.set_title("Pedersen")
# fig2.suptitle("Conductivities \nyear: 2008 lat: -25 lon: 125")

#====

fig3,(ax1,ax2) = plt.subplots(2,1,figsize=(10, 20))
# for i in range(len(hour)):
#     ax1.semilogx(condHall2000[i]*-1,cond2000[i].iri.dado["H(km)"],"--",label = hour[i]+" year 2000")
#     ax1.semilogx(condHall2008[i]*-1,cond2008[i].iri.dado["H(km)"],label = hour[i]+" year 2008")
# ax1.legend() 
# # ax1.set_title("Hall")

dif_condHall = [condHall2000[i]*(-1) - condHall2008[i]*(-1) for i in range(len(hour))]
# i=0
# ax1.plot(dif_condHall[i],cond2000[i].iri.dado["H(km)"],label = hour[i] + "2000-2008")
# ax1.plot(condHall2000[i],cond2000[i].iri.dado["H(km)"],label = hour[i]+" 2000")
# ax1.plot(condHall2008[i],cond2000[i].iri.dado["H(km)"],label = hour[i]+" 2008")
# ax1.grid()
# ax2.plot(condPedersen2000[i] - condPedersen2008[i],cond2008[i].iri.dado["H(km)"],"--",label = hour[i])
# ax1.legend() 

for i in range(len(hour)):
    ax1.plot(dif_condHall[i],cond2000[i].iri.dado["H(km)"],label = hour[i] + "2000-2008")
    # ax1.plot(condHall2000[i],cond2000[i].iri.dado["H(km)"],label = hour[i]+" 2000")
    # ax1.plot(condHall2008[i],cond2000[i].iri.dado["H(km)"],label = hour[i]+" 2008")

    ax2.plot(condPedersen2000[i] - condPedersen2008[i],cond2008[i].iri.dado["H(km)"],"--",label = hour[i]+"2000-2008")
#ax1.set_xscale('log')
ax2.legend() 
ax1.legend() 
ax1.grid()
# for i in range(len(hour)):
#     ax1.semilogx(dif_condHall[i],cond2000[i].iri.dado["H(km)"],label = hour[i])
#     ax2.semilogx(condPedersen2000[i] - condPedersen2008[i],cond2008[i].iri.dado["H(km)"],"--",label = hour[i])
ax1.legend() 
ax1.set_title("Hall 2000-2008")
ax2.set_title("Pedersen 2000-2008")
