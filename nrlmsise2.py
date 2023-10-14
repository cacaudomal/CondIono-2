#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:17:20 2022

@author: clara
"""

import pandas as pd
import matplotlib.pyplot as plt

class nrlmsise():
    def __init__(self,nome_arq):
        #self.dado = self.pegar_dado(nome_arq)
#        self.cabecalho = self.separar_cabecalho()
        self._read_data(nome_arq)
        self._convert_cm3tom3()
        #print("\n\n class = nrlmsise - __init__ : Cabecalho:\n",self.cabecalho)
        #print("\n class : nrlmsise - __init__  - Dado:\n",self.dado)
    
    
    def pegar_dado(self,nome_arq):
        return pd.read_fwf(nome_arq,
                   names = ["Year","Mon","Day","DOY","hour","H(km)","Lat","Lon","O (cm-3)","N2 (cm-3)","O2 (cm-3)","air(gm/cm3)","Tn(K)","exoT(K)","He (cm-3)","Ar (cm-3)","H (cm-3)","N (cm-3)","F107", "F107a","apdaily","ap0-3","ap3-6","ap6-9","ap9-12","ap12-33","ap33-59"])
    
    
    def _read_data(self,nomearq):
         """
         FUNÇÃO QUE LÊ O ARQUIVO DO IRI 2020, SALVA CADA LINHA NUMA LISTA, RETORNA 
         OS VALORES A PARTIR DA LINHA 32 E DESCARTA O RESTO.
    
         Parameters
         ----------
         nomearq : STRING
             NOME DO ARQUIVO ONDE ESTÃO OS DADOS. IRI2020
    
         Returns
         -------
         guardadado : LIST
             LISTA DE STRINGS CONTENDO OS DADOS EM CADA LINHA. \
                 PULAMOS PARA A LINHA 32 PQ É ONDE ESTÁ O NOME DAS COLUNAS. NA LINHA 33 OS DADOS COMEÇAM DE FATO.
    
         """
         with open(nomearq,"r") as arq: 
             dado = arq.readlines() #le todas as linhas do arquivo
        
         guarda_dado = []
         
         for i in range(len(dado)):
             guarda_dado.append(dado[i].split())#separa os valores no string da linha pro seus próprios espaços
         
         dado_float = [list(map(float,i)) for i in guarda_dado]
         
         dado = pd.DataFrame(dado_float, columns =  ["Year","Mon","Day","DOY","hour","H(km)","Lat","Lon","O (cm-3)","N2 (cm-3)","O2 (cm-3)","air(gm/cm3)","Tn(K)","exoT(K)","He (cm-3)","Ar (cm-3)","H (cm-3)","N (cm-3)","F107", "F107a","apdaily","ap0-3","ap3-6","ap6-9","ap9-12","ap12-33","ap33-59"])   
         
         self.dado = dado
        
    def _convert_cm3tom3(self):
        #j=("O (m^3)","N2 (m^3)","O2 (m^3)")
        for i in ("O (cm-3)","N2 (cm-3)","O2 (cm-3)","He (cm-3)","Ar (cm-3)","H (cm-3)","N (cm-3)"):
            #print("\n\nfora SI:\n",self.dado[i])
            self.dado[i] = self.dado[i] * 1e6
            #print("\nno SI",self.dado[i])
            
        self.dado.columns = ["Year","Mon","Day","DOY","hour","H(km)","Lat","Lon","O (m-3)","N2 (m-3)","O2 (m-3)","air(gm/cm3)","Tn(K)","exoT(K)","He (m-3)","Ar (m-3)","H (m-3)","N (m-3)","F107", "F107a","apdaily","ap0-3","ap3-6","ap6-9","ap9-12","ap12-33","ap33-59"]
    
    
    def separar_cabecalho(self):
        valores = []
        colunas_repetidas = []

        colunas = ["Year","Mon","Day","DOY","hour","H(km)","Lat","Lon","F107", "F107a","apdaily","ap0-3","ap3-6","ap6-9","ap9-12","ap12-33","ap33-59"] 
        
        for col in colunas:
            if self.dado[col].min() == self.dado[col].max(): #CHECAMOS AS COLUNAS QUE TEM O VALOR MAXIMO E MÍNIMO IGUAL
                valores.append(self.dado[col][1])
                colunas_repetidas.append(col)
        
        self.dado = self.dado.drop(columns = colunas_repetidas) #tirando as colunas repetidas da tabela original

        cabecalho = pd.Series(valores,index = colunas_repetidas)
        
        return cabecalho
    
    
    def plot_concentracoes(self):
        #plt.figure(figsize=(5,5))
        H = self.dado["H(km)"]
        lat = list(self.dado["Lat"])
        lon = list(self.dado["Lon"])
        hour = list(self.dado["hour"])
        fig, ax = plt.subplots(figsize=(10, 5))
        
        ax.semilogx(self.dado["O (m-3)"],H,'-', label='O')
        ax.semilogx(self.dado["O2 (m-3)"],H,'-', label='$O_2$')
        ax.plot(self.dado["N2 (m-3)"],H, label='$N_2$')
        ax.semilogx(self.dado["H (m-3)"],H,'-', label='H')
        ax.semilogx(self.dado["He (m-3)"],H,'-', label='He')
        #ax.semilogx(self.dado["N (cm-3)"],self.dado["Heigth(km)"],'-', label='N')
        ax.semilogx(self.dado["Ar (m-3)"],H,'-', label='Ar')

        ax.set_title('Perfil da Composição da Atmosfera Neutra\n lat: '+str(lat[0])+" lon:"+str(lon[0])+" hour: "+str(hour[0]))
        
        ax.set_xlabel("$log_{10}$ Densidade ($m^{-3}$)")
        ax.set_ylabel("Altura (km)")
         
        #plt.xlim(left=10e2)
        
        ax.legend()
        ax.grid(True)
        
        #ax.show()


# nomearqMSISE = "nrlmsise_lat69_lon19_z300_s5_data30MAR2012_h6UT_2.txt"
# msise = nrlmsise(nomearqMSISE)
# msise.plot_concentracoes()

# nomearqMSISE = "NRLMSISE2_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt"
# msise_b =  nrlmsise(nomearqMSISE)
# msise_b.plot_concentracoes()

#print(a.dado)