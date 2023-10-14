# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 19:26:31 2023

@author: tedea
"""
import matplotlib.pyplot as plt
import pandas as pd
import pyigrf_clara_0_4 as igrf
import condutividade_0_7 as cd

#====================== TESTE ====================================
name_saida = "teste_cond"
nomearqIRI12LT = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt"
nomearqMSISE12LT = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt"

nomearqIRI00LT = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt"
nomearqMSISE00LT = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt"

nomearqIRI18LT = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h18LT.txt"
nomearqMSISE18LT = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h18LT.txt"

nomearqIRI06LT = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h06LT.txt"
nomearqMSISE06LT = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h06LT.txt"

igrf_2008 = igrf.IGRF(-25,-55,0,2008,name_saida)
igrf_2008.calc_perfilh(-25,55,lim_h = 755)

igrf_2000 = igrf.IGRF(-25,55,0,2000,name_saida)
igrf_2000.calc_perfilh(-25,55,lim_h = 755)

#calculando as condutividades
cond_12LT = cd.condiono_adachi(igrf_2000.DfPerfilh["B(T)"], nomearqIRI12LT, nomearqMSISE12LT)
cond_00LT = cd.condiono_adachi(igrf_2000.DfPerfilh["B(T)"], nomearqIRI00LT, nomearqMSISE00LT)
cond_18LT = cd.condiono_adachi(igrf_2000.DfPerfilh["B(T)"], nomearqIRI18LT, nomearqMSISE18LT)
cond_06LT = cd.condiono_adachi(igrf_2000.DfPerfilh["B(T)"], nomearqIRI06LT, nomearqMSISE06LT)

cond_12LT.calc_Hall(cond_12LT.Freq["fen"],cond_12LT.Freq["fin1"],cond_12LT.Freq["fin2"],\
                    cond_12LT.girofreq["wi1"],cond_12LT.girofreq["wi2"],cond_12LT.girofreq["we"],\
                    cond_12LT.p1, cond_12LT.p2, cond_12LT.iri.dado["Ne(m-3)"], cond_12LT.B)
cond_12LT.calc_Pedersen(cond_12LT.Freq["fen"],cond_12LT.Freq["fin1"],cond_12LT.Freq["fin2"],\
                 cond_12LT.girofreq["wi1"],cond_12LT.girofreq["wi2"],cond_12LT.girofreq["we"],\
                 cond_12LT.p1, cond_12LT.p2, cond_12LT.iri.dado["Ne(m-3)"], cond_12LT.B)

cond_00LT.calc_Hall(cond_00LT.Freq["fen"],cond_00LT.Freq["fin1"],cond_00LT.Freq["fin2"],\
                    cond_00LT.girofreq["wi1"],cond_00LT.girofreq["wi2"],cond_00LT.girofreq["we"],\
                    cond_00LT.p1, cond_00LT.p2, cond_00LT.iri.dado["Ne(m-3)"], cond_00LT.B)
cond_00LT.calc_Pedersen(cond_00LT.Freq["fen"],cond_00LT.Freq["fin1"],cond_00LT.Freq["fin2"],\
                 cond_00LT.girofreq["wi1"],cond_00LT.girofreq["wi2"],cond_00LT.girofreq["we"],\
                 cond_00LT.p1, cond_00LT.p2, cond_00LT.iri.dado["Ne(m-3)"], cond_00LT.B)
    

cond_18LT.calc_Hall(cond_18LT.Freq["fen"],cond_18LT.Freq["fin1"],cond_18LT.Freq["fin2"],\
                    cond_18LT.girofreq["wi1"],cond_18LT.girofreq["wi2"],cond_18LT.girofreq["we"],\
                    cond_18LT.p1, cond_18LT.p2, cond_18LT.iri.dado["Ne(m-3)"], cond_18LT.B)
cond_18LT.calc_Pedersen(cond_18LT.Freq["fen"],cond_18LT.Freq["fin1"],cond_18LT.Freq["fin2"],\
                 cond_18LT.girofreq["wi1"],cond_18LT.girofreq["wi2"],cond_18LT.girofreq["we"],\
                 cond_18LT.p1, cond_18LT.p2, cond_18LT.iri.dado["Ne(m-3)"], cond_18LT.B)

cond_06LT.calc_Hall(cond_06LT.Freq["fen"],cond_06LT.Freq["fin1"],cond_06LT.Freq["fin2"],\
                    cond_06LT.girofreq["wi1"],cond_06LT.girofreq["wi2"],cond_06LT.girofreq["we"],\
                    cond_06LT.p1, cond_06LT.p2, cond_06LT.iri.dado["Ne(m-3)"], cond_06LT.B)
cond_06LT.calc_Pedersen(cond_06LT.Freq["fen"],cond_06LT.Freq["fin1"],cond_06LT.Freq["fin2"],\
                 cond_06LT.girofreq["wi1"],cond_06LT.girofreq["wi2"],cond_06LT.girofreq["we"],\
                 cond_06LT.p1, cond_06LT.p2, cond_06LT.iri.dado["Ne(m-3)"], cond_06LT.B)
    

CondH12LT = [cond_12LT.CondH[i] * -1 for i in range(len(cond_12LT.CondH))]
CondH00LT = [cond_00LT.CondH[i] * -1 for i in range(len(cond_00LT.CondH))]

CondH18LT = [cond_12LT.CondH[i] * -1 for i in range(len(cond_12LT.CondH))]
CondH06LT = [cond_06LT.CondH[i] * -1 for i in range(len(cond_06LT.CondH))]

#plotando
fig,(ax,ax1,ax2) = plt.subplots(1,3,figsize=(15, 5))

ax.semilogx(CondH00LT[17:],range(80,750,5),label="00h LT")
ax.semilogx(CondH06LT[17:],range(80,750,5),label="06h LT")
ax.semilogx(CondH12LT[17:],range(80,750,5),label="12h LT")
ax.semilogx(CondH18LT[17:],range(80,750,5),label="18h LT")
ax.grid()
ax.set_title("Hall")
ax.legend()


ax1.semilogx(cond_00LT.CondP,range(0,755,5),label="00h LT")
ax1.semilogx(cond_06LT.CondP,range(0,755,5),label="06h LT")
ax1.semilogx(cond_12LT.CondP,range(0,755,5),label="12h LT")
ax1.semilogx(cond_18LT.CondP,range(0,755,5),label="18h LT")

ax1.grid()
ax1.legend()
ax1.set_title("Pedersen")

ax2.semilogx(cond_00LT.iri.dado["Ne(m-3)"],cond_00LT.iri.dado["H(km)"],label="00LT")
ax2.semilogx(cond_06LT.iri.dado["Ne(m-3)"],cond_06LT.iri.dado["H(km)"],label="06LT")
ax2.semilogx(cond_12LT.iri.dado["Ne(m-3)"],cond_12LT.iri.dado["H(km)"],label="12LT")
ax2.semilogx(cond_18LT.iri.dado["Ne(m-3)"],cond_18LT.iri.dado["H(km)"],label="18LT")
ax2.legend()

fig.suptitle("Conductivities\n Year: 2000 lat: -25 lon: 125")


