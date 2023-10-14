# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:07:36 2021

@author: Clara

contém as classes freqcol e fen
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pyigrf_clara_0_4 as igrf
import nrlmsise2 as msise2
import iri2_0_2 as iri

class freqcol():
    """
    Classe para calcular as frequências de colisão dadas segundo Adachi et. al. 2017.
    
    ...
    
    Dependencies: 
    -------------
        : numpy, matplotlib, pandas       
        
    Attributes
    ----------
    nN2 : LIST FLOAT
        densidade de N2, [m^3] NRLMSISE2
    nO2 : LIST FLOAT
        densidade de O2 [m^3] NRLMSISE2
    nO : LIST FLOAT
        densidade de O [m^3] NRLMSISE2    
    Ti : PANDA SERIES
        TEMPERATURA DOS ÍONS EM [K].
    Tn : PANDA SERIES
        TEMPERATURA DAS PARTÍCULAS NEUTRAS EM [K].
    h : PANDA SERIES
        ALTURAS NAS QUAIS OS VALORES SERÃO CALCULADOS [km].
        
    Methods
    ----------
        set_Temp(self,Te,Tn,Ti)
        set_Tr(self,Ti,Tn)
        set_atmNeutra(self,nN2,nO2,nO)
        calc_freq(self,h)
        calc_fen(self)
        plot_freq_h(self)
        
    """
    
    def __init__(self, nN2, nO2, nO,Te,Tn,Ti,h):
        self.set_atmNeutra(nN2, nO2, nO)
        self.set_temp(Te,Tn,Ti)
        self.set_Tr(Ti,Tn)
        
        self.calc_freq(h)
        
    def set_temp(self,Te,Tn,Ti):
        """
        SETS THE TEMPERATURES VALUES FOR THE ATMOSPHERE.

        Parameters
        ----------
        Te : PANDA SERIES
            TEMPERATURA DOS ELÉTRONS EM [K].
        Ti : PANDA SERIES
            TEMPERATURA DOS ÍONS EM [K].
        Tn : PANDA SERIES
            TEMPERATURA DAS PARTÍCULAS NEUTRAS EM [K].

        Returns
        -------
        None.

        """
        self.Te = Te
        self.Ti = Ti    
        self.Tn = Tn
    
    
    def set_Tr(self,Ti,Tn):
        """
        CRIA UMA PANDA SERIE COM OS VALORES DA TEMPERATURA RELATIVA E PREENCHE SEUS VALORES.

        Parameters
        ----------
        Ti : PANDA SERIES
            TEMPERATURA DOS ÍONS EM [K].
        Tn : PANDA SERIES
            TEMPERATURA DAS PARTÍCULAS NEUTRAS EM [K].

        Returns
        -------
        None.

        """
        self.Tr = (Ti + Tn)/2  
        
        
    def set_atmNeutra(self,nN2,nO2,nO):
        """
        SETS THE VALUES FOR nN2, nO2 and nO.

        Parameters
        ----------
        nN2 : PANDA SERIES - FLOATS
            DENSITY OF N2 AT A GIVEN HEIGHT [m^3] NRLMSISE2
        nO2 : PANDA SERIES - FLOATS
            DENSITY OF O2 AT A GIVEN HEIGHT [m^3] NRLMSISE2
        nO : PANDA SERIES - FLOATS
            DENSITY OF O AT A GIVEN HEIGHT [m^3] NRLMSISE2

        Returns
        -------
        None.

        """
        self.nN2 = nN2
        self.nO2 = nO2 #lista com a densidade de O2 [cm^-3] é convertido para [m^-3]
        self.nO  = nO #densidade de O [cm^-3] é convertido para [m^-3]
        
        
    def calc_freq(self,h):
        """
        Cálculo das frequêncisa de colisão dos íon O+, ion 1 e elétrons com as partículas neutras.        
        Equações retiradas do trabalho de Adachi et. al, 2017.
    
        Parameters
        ----------
        h : DATAFRAME
            ALTURA PARA QUAL A FREQUÊNCIA DE COLISÃO SERÁ CALCULADA [km]
            
        Returns
        -------
            freqcol : DATAFRAME 
                UM DATAFRAME CUJAS COLUNAS SÃO RESPECTIVAMENTE:
                fen : frequência de colisão dos elétrons com as partículas neutras [Hz]
                fin1 : frequência de colisão do íon 1 com as partículas neutras [Hz]
                fin2 : frequência de colisão do íon 2 com as partículas neutras [Hz]
                H(km) : ALTURA PARA QUAL A FREQUÊNCIA DE COLISÃO FOI CALCULADA [km]
        """
        
        fen = 2.33e-17 * self.nN2 * (1 - 1.21e-4*self.Te) * self.Te + 1.82e-16 * self.nO2 * (1 + 3.6e-2*np.sqrt(self.Te)) * np.sqrt(self.Te) + 8.9e-17*self.nO * (1 + 5.7e-4*self.Te) * np.sqrt(self.Te)
        fin1 = (4.29 * self.nN2 + 4.23 * self.nO2 + 2.41 * self.nO) * 1e-16
        fin2 = 6.82e-16*self.nN2 + 6.66e-16*self.nO2 + 3.32e-17*self.nO * np.sqrt(self.Tr) * (1.08 - 0.139*np.log10(self.Tr) + 4.51e-3*(np.log10(self.Tr)**2))
        #print("\n\ncalc_freq- fen:",type(fen))
        freqcol =  pd.concat([fen,fin1,fin2,h], axis=1, keys=["fen","fin1","fin2","H(km)"])
        
        self.freqcol = freqcol
        #print("freqcol - calc_freq type freqcol:")
        #print("\nFREQCOL - calc_freq- freqcol: ",self.freqcol)
        
        return freqcol
    
    
    def plot_freq_h(self):
        """
        CRIA PERFIL DE ALTURA DAS FREQUÊNCIAS DE COLISÕES.

        Returns
        -------
        None.

        """
        fig,ax = plt.subplots()
        
        ax.semilogx(self.freqcol["fen"],self.freqcol["H(km)"],'-', label='fen (Hz)')
        ax.semilogx(self.freqcol["fin1"],self.freqcol["H(km)"],'-', label='fin1 (Hz)')
        ax.semilogx(self.freqcol["fin2"],self.freqcol["H(km)"],'-', label='fin2 (Hz)')
        
        ax.set_title("Frequências de Colisão")
        ax.set_xlabel("$log_{10}$ Frequência de Colisão ($Hz$)")
        ax.set_ylabel("Altura (km)")
        
        ax.legend()
        ax.grid(True)
    
    
#======= TESTE ==========================
# nomearqMSISE = "nrlmsise_lat69_lon19_z300_s5_data30MAR2012_h6UT.txt"
# nomearqIRI = "IRI2016_lat69_lon19_z300_s5_data30MAR2012_h6UT.txt"

# tmsise = msise2.nrlmsise(nomearqMSISE)
# tiri = iri.iri(nomearqIRI)

# freq = freqcol(tmsise.dado["N2 (m-3)"], tmsise.dado["O2 (m-3)"], tmsise.dado["O (m-3)"], tiri.dado["Te/K"], tiri.dado["Tn/K"], tiri.dado["Ti/K"],tiri.dado["H(km)"])

