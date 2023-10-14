# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 12:19:46 2021

@author: Clara
"""

import iri2_0_1 as iri
import pyigrf_clara_0_3 as igrf
import nrlmsise2 as msise
import freqcol_0_3_1 as fc
import condutividade_0_5 as cond
import matplotlib.pyplot as plt

#import numpy as np
saida_igrf = "SAIDA_2008_final.txt"
nomearqIRI = "IRI_lat-25_lon125_z0_750_s_5data01JAN2008_h00LT.txt"
nomearqMSISE = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2008_h00LT.txt"

#a = d.dado("dadosIGRF.txt")

v_iri = iri.iri(nomearqIRI)

v_igrf = igrf.IGRF(0.0,50,100,2008,saida_igrf)
#v_igrf.calc_grid(intervalo_h = 10,lim_h = 300,intervalo_lat=-10,lim_lat=-60,intervalo_lon=-10,lim_lon=-120)
v_igrf.calc_perfilh(-25,-55)
v_msise = msise.nrlmsise(nomearqMSISE)
v_freqcol = fc.freqcol(v_msise.dado['N2 (cm-3)'],v_msise.dado['O2 (cm-3)'],v_msise.dado['O (cm-3)'],v_iri.dado['Te/K'],v_iri.dado['Tn/K'],v_iri.dado['Ti/K'], v_iri.dado["H(km)"])

#centro_anom = v_igrf.coord_centro(v_igrf.Dfgrid)

condutividade = cond.condiono_adachi(v_igrf.DfPerfilh["B(T)"], nomearqIRI, nomearqMSISE)
condutividade.plot_Hall( v_iri.dado["H(km)"])
#v_igrf.calc_perfilh(500,10)
#v_igrf.calc_grid(intervalo_h = 10,lim_h = 300,intervalo_lat=-10,lim_lat=-60,intervalo_lon=-10,lim_lon=-120)
#plt.scatter(v_igrf.DfPerfilh['Total_intensity'],v_igrf.DfPerfilh['Altitude'])
#plt.xlabel("Total Intensity")
#plt.ylabel("H (km)")
#print("c:\n",v_igrf.parametros)

#plt.plot(v_iri.dado["H(km)"],v_fen.fen)
#plt.plot(v_fin1.fin1,v_iri.dado["H(km)"])
#plt.plot(v_fin2.fin2,v_iri.dado["H(km)"])
# with open("arq-test1.txt","w") as file:
#     i=0
#     while i < len(g.fen):
#         file.write(str(g.fen[i])+"\n")
#         i+=1
        
#print("d.__doc__:",g.__doc___)