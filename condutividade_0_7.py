# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 11:42:46 2023

@author: tedea
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


import pyigrf_clara_0_4 as igrf
import freqcol_0_3_1 as fc
import nrlmsise2 as msise2
import iri2_0_2 as iri

class condiono_adachi():
    def __init__(self,B,nomearqIRI,nomearqNRLMSISE2):
        self.me = 9.109389e-31 #Massa do elétron em repouso [kg]
        self.e = -1.602177e-19 #Carga do elétron [C]
        self.mi1 = 5.065e-26 #Massa do íon 1 uma mistura de NO+ (75%) e O2+ (25%) (30.5 u.m.a.) [kg]
        self.mi2 = 2.657e-26 #Massa do íon 2 (O+) [kg] (16 a.m.u)
        
        #self._set_massaIon()
        self._set_nomearq(nomearqIRI, nomearqNRLMSISE2) #salvando o nome dos arquivos em variaveis na classe
        
        self._set_iri(nomearqIRI)
        self._set_nrlmsise2(nomearqNRLMSISE2)
        self._set_B(B)
        
        self.calc_freqcol(self.msise2.dado['N2 (m-3)'],self.msise2.dado['O2 (m-3)'],\
                          self.msise2.dado['O (m-3)'],self.iri.dado['Te/K'],\
                          self.iri.dado['Ti/K'],self.iri.dado['Tn/K'],self.iri.dado["H(km)"])
        
        self.calc_all_girofreq(self.B)
        self.calc_prelativa_all(self.iri.rhodado["O+"], self.iri.rhodado["NO+"],\
                                self.iri.rhodado["O2+"],self.iri.dado["Ne(m-3)"])
        #print("self.Freq[fen]",self.Freq["fen"]) 
        
        
    def _set_massaIon(self):
        self.mi1 = 5.065e-26 #é a massa do íon 1 (30.5 u.m.a.) em kg
        self.mi2 = 2.657e-26 #é a massa do íon 2 (O+) em kg
        
        
    def _set_nomearq(self,nomearqIRI,nomearqNRLMSISE2):
         self.nomearqIRI = nomearqIRI
         self.nomearqMSISE = nomearqNRLMSISE2
     
        
    def _set_iri(self,nomearqIRI):   
        self.iri = iri.iri(nomearqIRI)
     
        
    def _set_nrlmsise2(self,nomearqMSISE):
        self.msise2 = msise2.nrlmsise(nomearqMSISE)
        #print("\n\n_set_nrlmsise2 - msise2_dado:",self.msise2.dado)
        
        
    def _set_B(self,B):
        """
        SETS THE TOTAL INTENSITY OF THE MAGNETIC FIELD.

        Parameters
        ----------
        B : PANDA SERIES
            TOTAL INTENSITY OF THE MAGNETIC FIELD [T].

        Returns
        -------
        None.

        """
        self.B = B     


    def calc_girofreq(self,mi,B):
        """
        Funcao para cálculo da girofrequência ou frequência de ciclotron da partícula num ponto.
        
        Parameters
        ----------
        B : FLOAT
            campo magnético em [T]
        mi : FLOAT
            massa do íon ou do elétron [kg]    
        q : FLOAT
            carga do íon/elétron [C]
        
        Returns:
        ----------
        wi : FLOAT
            girofrequencia [Hz]
            
        """
        wi = np.sqrt(self.e**2) * np.sqrt(B**2)/mi
        
        return wi 
    
    
    def calc_all_girofreq(self,B):
        print("calculando as freqcol all: start ")
        wi1 = self.calc_girofreq(self.mi1, B)
        wi2 = self.calc_girofreq(self.mi2, B)
        we = self.calc_girofreq(self.me, B)
        #print("\n\ncondiono_adachi calc_all_girofreq  wi1: ",wi1,"\nwi2:",wi2,"\nwe5:",we)
        
        self.girofreq = pd.concat([we,wi1,wi2], axis=1, keys=["we","wi1","wi2"])
           
        
    def _calc_pRelativa(self,rhoi,ne):
        """
        Calcula a densidade numérica relativa da espécie ionica. Brekke (1993)
        
        Parameters
        ----------
        rhoi : PANDA SERIES FLOAT
            densidade do íon [m^-3]     
        ne : PANDA SERIES FLOAT
            Densidade de elétrons [elétrons/m^3]
        
        Returns:
        ----------
        pi : PANDA SERIES FLOAT
            densidade numérica relativa (Brekke,1983)
            
        """    
        pi = rhoi/ne        
        return pi
    
    
    def calc_prelativa_all(self,rho_íonO,rho_íonNO,rho_íonO2,ne):
        self.rho1,self.rho2 = self.calc_rho_numion(rho_íonO,rho_íonNO,rho_íonO2,ne)
        
        self.p1 = self._calc_pRelativa(self.rho1, ne)
        self.p2 = self._calc_pRelativa(self.rho2, ne)
        
        return self.p1,self.p2
           
    
    def calc_rho_numion(self,rho_íonO,rho_íonNO,rho_íonO2,ne):
        """
        Calcula a Densidade numérica dos íons O+ e fictício 1 (Brekke,1983).
        Razão entre o número de íons e o volume.
        
        Parameters
        ----------
        rho_íonO : LIST FLOAT
            concentração do íon O+ [%]

        rho_íonNO : LIST FLOAT
            concentração do íon NO+ [%]

        rho_íonO2 : LIST FLOAT
            concentração do íon O2+ [%]

        ne : LIST FLOAT
            Densidade de elétrons [elétrons/m^3]
        
        Returns:
        ----------
        rhoi1  : FLOAT
            densidade do íon fictício 1 [m^-3]
        rhoi2 : FLOAT
            densidade do íon O+ [m^-3]    
        """
        
        #densidade do íon ficticio 1 [m^-3]
        rhoi1 = ne * (rho_íonNO + rho_íonO2)/100
        #print(self.h[i],"rho1",rhoi1,"m^-3")
    
        #densidade do íon O+ [m^-3]
        rhoi2 = ne * rho_íonO/100 
        #print("rhoi2",rhoi2,"m^-3\n")
        return rhoi1, rhoi2
    
    
    def calc_freqcol(self,rhoN2,rhoO2,rhoO,Te,Ti,Tn,h):
        '''
        
        CALCULA AS FREQUÊNCIAS DE COLISÃO.

        Parameters
        ----------
        rhoN2 : PANDA SERIES - FLOAT
           DENSITY OF N2 AT A GIVEN HEIGHT [m^3] NRLMSISE2
        rhoO2 : PANDA SERIES - FLOAT
            DENSITY OF O2 AT A GIVEN HEIGHT [m^3] NRLMSISE2
        rhoO : PANDA SERIES - FLOAT
            DENSITY OF O AT A GIVEN HEIGHT [m^3] NRLMSISE2
            
        Te : PANDA SERIES - FLOAT
            TEMPERATURA DOS ELÉTRONS [K].
        Ti : PANDA SERIES - FLOAT
            TEMPERATURA DOS ÍONS [K].
        Tn : PANDA SERIES - FLOAT
            TEMPERATURA DAS PARTÍCULAS NEUTRAS [K].

        Returns
        -------
        
        '''
        a = fc.freqcol(rhoN2, rhoO2, rhoO, Te, Tn, Ti,h)
        self.Freq = a.calc_freq(h)
        
        #print("\n\Condutividade_0_7 - type freq : ",type(self.Freq),"\n")
        return self.Freq
    
    
    def calc_Hall(self,fen,fin1,fin2,wi1,wi2,we,p1,p2,ne,B):
        """
        CALCULA A CONDUTIVIDADE DE HALL APARTIR DAS EQUAÇÕES DE Adachi et al.
        Earth, Planets and Space (2017).

        Parameters
        ----------
        fen : PANDA SERIES
            frequência de colisão dos elétrons com as partículas neutras [Hz].
        fin1 : PANDA SERIES
            frequência de colisão do íon 1 com as partículas neutras [Hz].
        fin2 : PANDA SERIES
            frequência de colisão do íon 2 com as partículas neutras [Hz].
        wi1 : PANDA SERIES
            girofrequência do íon 1 [Hz].
        wi2 : PANDA SERIES
            girofrequência do íon 2 [Hz].
        we : PANDA SERIES
            girofrequência do elétron [Hz].
        p1 : TYPE
            DESCRIPTION.
        p2 : TYPE
            DESCRIPTION.
        ne : TYPE
            densidade de elétrons em [m^-3].
        B : TYPE
            intensidade do campo magnético da Terra [T].

        Returns
        -------
        cond_hall : TYPE
            DESCRIPTION.

        """         
        print("Calculando a Condutividade de Hall")
            
        a1 = (wi2**2)/(wi2**2 + fin2**2)
        b1 = (wi1**2)/(wi1**2 + fin1**2)
        c1 = (we**2)/(we**2 + fen**2)

        soma = c1 - (p1*b1) - (p2*a1)
        d = (ne * self.e)/B

        self.CondH  = d * soma
        
        return self.CondH   
       
        
    def calc_Pedersen(self,fen,fin1,fin2,wi1,wi2,we,p1,p2,ne,B):
        """
        # CALCULA A CONDUTIVIDADE DE PEDERSEN APARTIR DAS EQUAÇÕES DE Adachi et al.
        # Earth, Planets and Space (2017).

        # Parameters
        # ----------
 
        # fen : PANDA SERIES
        #     frequência de colisão dos elétrons com as partículas neutras [Hz].
        # fin1 : PANDA SERIES
        #     frequência de colisão do íon 1 com as partículas neutras [Hz].
        # fin2 : PANDA SERIES
        #     frequência de colisão do íon O+ com as partículas neutras [Hz].
        # wi1 : PANDA SERIES
        #     girofrequencia do íon 1 [Hz].
        # wi2 : TYPE
        #     girofrequencia do íon O+ [Hz].
        # we : TYPE
        #     DESCRIPTION.
        # p1 : TYPE
        #     DESCRIPTION.
        # p2 : TYPE
        #     DESCRIPTION.
        # ne : TYPE
        #     DESCRIPTION.
        # B : TYPE
        #     DESCRIPTION.
        # rho_íonNO : TYPE
        #     DESCRIPTION.

        # Returns
        # -------
        # condP : TYPE
        #     DESCRIPTION.

        # """
        d = 0
        soma = 0
        
        print("===== Calculando a condutividade de Pedersen =====")
        
       
        a1 = (wi2 * fin2)/(wi2**2 + fin2**2)
        b1 = (wi1 * fin1)/(wi1**2 + fin1**2)
        c1 = (we * fen)/(we**2 + fen**2)
            
        soma = c1 + p1 * b1 + p2 * a1
        d = (ne * np.sqrt(self.e**2))/B
        
        self.CondP = d * soma
        
        return self.CondP
    
    
    def plot_Hall(self):
        h = self.msise2.dado["H(km)"]
        plt.figure(figsize=(5,5))
        
        condH = [self.CondH[i] * -1 for i in range(len(self.CondH))]
        
        plt.plot(condH,h)
        plt.xscale('log')
        plt.ylabel("Height (km)")
        plt.xlabel("$log_{10}$ Conductivity (S/m)",fontsize=15)
        #ax.legend(title="Hour")
        plt.title("Hall conductivity")
        plt.grid()  
        
        plt.show()
          
        
    def plot_Pedersen(self):
        plt.figure(figsize=(5,5))
        
        plt.plot(self.CondP,self.msise2.dado["H(km)"],label="Pedersen Conductivity")
       
        plt.xscale('log')

        plt.ylabel("Height (km)")
        plt.xlabel("$log_{10}$ Conductivity (S/m)",fontsize=15)
        plt.title("Pedersen conductivity")  
        
        plt.grid()  
        
        plt.show()
       
        
    def plot_girofreq(self):
        
        h = self.msise2.dado["H(km)"]
        
        plt.figure(figsize=(5,5))

        #plt.plot(self.girofreq["we"],h,label="we")
        plt.plot(self.girofreq["wi1"],h,label="wi1")
        plt.plot(self.girofreq["wi2"],h,label="wi2")

        plt.title("girofrequencias com a altura (km)")
        plt.ylabel("Height (km)")
        plt.xlabel("$log_{10}$ Frequência de Ciclotron (Hz)")
        plt.legend()
        #plt.xscale('log')  

        plt.grid()
        plt.show()


    def plot_freqcol(self):
        h = self.Freq["H(km)"]
        
        plt.figure(figsize=(5,5))

        plt.plot(self.Freq["fen"],h,label="electron")
        plt.plot(self.Freq["fin1"],h,label="ficticious ion 1")
        plt.plot(self.Freq["fin2"],h,label="O+")
        
        plt.title('Collision Frequency')
        plt.xlabel("$log_{10}$ Collision Frequency (Hz)")
        plt.ylabel("Height (km)")
        
        plt.xscale('log')  

        plt.legend()
        plt.grid()
        
        plt.show()
        
    
#====================== TESTE ====================================
# name_saida = "grid_igrf_2000"
# nomearqIRI = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt"
# nomearqMSISE = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt"

# # nomearqMSISE = "nrlmsise_lat69_lon19_z300_s5_data30MAR2012_h6UT.txt"
# # nomearqIRI = "IRI2016_lat69_lon19_z300_s5_data30MAR2012_h6UT.txt"
# nomearqIRI2 = "IRI2016_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt"
# nomearqMSISE2 = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt"
# path = Path("Dado_IRI/IRI2020_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt")
# assert path.exists()
# #igrf_a = igrf.IGRF(-25,-55,0,2008,name_saida)
# #igrf_a.get_grid("teste_cond_grid.csv")

# igrf_b = igrf.IGRF(-25,55,0,2000,name_saida)
# igrf_b.calc_perfilh(-25,55,lim_h = 755)

# a = condiono_adachi(igrf_b.DfPerfilh["B(T)"], path, nomearqMSISE)

#igrf_a.calc_grid()


#centro_anom = igrf_a.coord_centro(igrf_a.Dfgrid)
#igrf_a.Dfgrid["Latitude"]

#d=(igrf_a.Dfgrid['Latitude'] ==	centro_anom[0])#pegando a latitude determinada pela função coord_centro

#=== Filtrando para pegar os dados na latitude e longitude do centro da anomalia =====

#c = igrf_a.Dfgrid.loc[(igrf_a.Dfgrid['Latitude'].round(4) == -25.0)]
#b_centro = igrf_a.Dfgrid.loc[(igrf_a.Dfgrid['Latitude'].round(4) == -25.0) & (igrf_a.Dfgrid['Longitude'].round(4) == -55.0)]
#b_centro = igrf_a.Dfgrid.loc[(igrf_a.DfPerfilh['Latitude'].round(4) == -25.0) & (igrf_a.DfPerfilh['Longitude'].round(4) == -55.0)]
#print("n\nblat25",b_lat_25_lon125)

# a = condiono_adachi(igrf_b.DfPerfilh["B(T)"], nomearqIRI, nomearqMSISE)
# a2 = condiono_adachi(igrf_b.DfPerfilh["B(T)"], nomearqIRI2, nomearqMSISE2)
# #a.msise2.plot_concentracoes()
# a.calc_Hall(a.Freq["fen"],a.Freq["fin1"],a.Freq["fin2"],\
#                     a.girofreq["wi1"],a.girofreq["wi2"],a.girofreq["we"],\
#                     a.p1, a.p2, a.iri.dado["Ne(m-3)"], a.B)
# #a.plot_Hall()

# a.calc_Pedersen(a.Freq["fen"],a.Freq["fin1"],a.Freq["fin2"],\
#                  a.girofreq["wi1"],a.girofreq["wi2"],a.girofreq["we"],\
#                  a.p1, a.p2, a.iri.dado["Ne(m-3)"], a.B)
# #a.plot_Pedersen()


# a2.calc_Hall(a2.Freq["fen"],a2.Freq["fin1"],a2.Freq["fin2"],\
#                     a2.girofreq["wi1"],a2.girofreq["wi2"],a2.girofreq["we"],\
#                     a2.p1, a2.p2, a2.iri.dado["Ne(m-3)"], a2.B)
# #a2.plot_Hall()

# a2.calc_Pedersen(a2.Freq["fen"],a2.Freq["fin1"],a2.Freq["fin2"],\
#                  a2.girofreq["wi1"],a2.girofreq["wi2"],a2.girofreq["we"],\
#                  a2.p1, a2.p2, a2.iri.dado["Ne(m-3)"], a2.B)
# #a2.plot_Pedersen()
# CondH1 = [a.CondH[i] * -1 for i in range(len(a.CondH))]
# CondH2 = [a2.CondH[i] * -1 for i in range(len(a2.CondH))]

# fig,(ax,ax1) = plt.subplots(1,2,figsize=(10, 5))
# ax.semilogx(CondH1[17:],range(80,750,5),label="12h LT")
# ax.semilogx(CondH2[17:],range(80,750,5),label="00h LT")
# ax.grid()
# ax.set_title("Hall")
# ax.legend()

# ax1.semilogx(a.CondP,range(0,755,5),label="12h LT")
# ax1.semilogx(a2.CondP,range(0,755,5),label="00h LT")
# ax1.grid()
# ax1.legend()
# ax1.set_title("Pedersen")
# fig.suptitle("Conductivities year: 2000 lat: -25 lon: 125")