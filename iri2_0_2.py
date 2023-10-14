#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 16:03:34 2022

@author: clara
"""
import matplotlib.pyplot as plt
import pandas as pd

class iri():
    """
    Classe para ler os dados do IRI2.
    
    ...
    
        
    Attributes
    ----------
    nome_arq : STRING
        NOME DO ARQUIVO IRI ONDE ESTÃO GUARDADOS OS DADOS A SEREM LIDOS
    
    Methods
    ----------
        cria_dado_dataframe(self)
        _ler_arq(self,nomearq)
        _string_para_float(self,dado)
    """
    
    def __init__(self,nome_arq):
        self.nome_arq = nome_arq
        
        #arqlido = self._ler_arq(self.nome_arq)
        #self.splitdata = self._string_para_float(arqlido[1:])
        self.cria_dado_dataframe()
        
        self.to_percentage()
        self.eletron_cm3to_m3()
        #self.rhonumion_all(self.dado["Ne(m-3)"])
        #self.cabecalho = self._pega_cabecalho()
        
        print("IRI - dado lido com sucesso.")
        
    
    def _ler_arq(self,nomearq):
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
            dado = arq.readlines() #le todas as linas do arquivo
        #print("ler_arq dado:\n",dado)
        
        guardadado = []
        for i in range(len(dado)):
            guardadado.append(dado[i].split())#separa os valores no string da linha pro seus próprios espaços
        #print("\nle iri : ",guardadado[32])
        
        return guardadado[31:] # 32 para IRI2020 e 31 para IRI2016
                               #33 é a ultima linha do cabeçalho a linha 32 tem o header do dataframe 
    
    def _string_para_float(self,dado):
        """
        Método para separar uma lista de strings contendo só números em uma lista de listas de floats.
        
        Parameters
        ----------
        dado : VETOR DE STRINGS
            DADOS A SEREM CONVERTIDOS.

        Returns
        -------
        result : VETOR DE FLOAT
            Dado convertido para float.  
            
        """
        result = [list(map(float,i)) for i in dado] #faz o typecast de todos os valores de cada linha da lista de listas de str para float
        
        return result    
         
    
    def cria_dado_dataframe(self):
        arqlido = self._ler_arq(self.nome_arq)
        splitdata = self._string_para_float(arqlido[1:])
        
        self.dado = pd.DataFrame(splitdata, columns = ["H(km)","Ne(cm-3)","Ne/NmF2","Tn/K","Ti/K","Te/K","O+","N+","H+","He+","O2+","NO+","Clust","TEC","t/%"])
        
    
    def calc_rho_numion(self, ne, rho_ion):
       """
       Calcula a Densidade numérica do ion a partir da densidade de elétrons e da concentração do íon na atmosfera.
    
       Parameters
       ----------
       rho_ion : PANDA SERIES - FLOAT
           concentração do íon [%]
       ne : PANDA SERIES - FLOAT
           Densidade de elétrons [elétrons/m^3]
    
       Returns:
       ----------
       rhoion : PANDA SERIES - FLOAT
           densidade numérica do íon [m^-3]
       
       """
       rhoion =  ne * rho_ion/100 #densidade do íon em [m^-3]
     
       return rhoion
 
    
    def eletron_cm3to_m3(self):
        """
        CRIA COLUNA COM OS VALORES DA DENSIDADE DE ELÉTRONS em M-³ NO DATAFRAME DADO.

        Returns
        -------
        None.

        """
        self.dado["Ne(m-3)"] = self.dado["Ne(cm-3)"] * 1e6
        
        
    def to_percentage(self):
        '''
        CONVERTE A CONCENTRAÇÃO DO ION DE PERCENTAGEM * 10 PARA PORCENTAGEM E COLOCA NUM DATAFRAME
        E SALVA NUM DATAFRAME SEPARADO.
        Returns
        -------
        None.

        '''
        self.rhodado = self.dado[["O+","N+","H+","He+","O2+","NO+"]] * 0.1


    def rhonumion_all(self,Ne):
        
        #Ne = self.dado["Ne(m-3)"]
        #x = pd.DataFrame(columns = ["O+","N+","H+","He+","O2+","NO+"])
        a = self.calc_rho_numion(Ne,self.rhodado["O+"])
        b = self.calc_rho_numion(Ne,self.rhodado["N+"])
        c = self.calc_rho_numion(Ne,self.rhodado["H+"])
        d = self.calc_rho_numion(Ne,self.rhodado["He+"])
        e = self.calc_rho_numion(Ne,self.rhodado["O2+"])
        f = self.calc_rho_numion(Ne,self.rhodado["NO+"])
        #print("\n\n\nabcd:",a,b,c,d,e,f)
        self.rhonum_dado = pd.concat([a,b,c,d,e,f], axis=1, keys=["O+","N+","H+","He+","O2+","NO+"])
        
        return self.rhodado
        
    
    def plot_densidade_e(self, ne, h, data=""):
        """
        FUNÇÃO PARA PLOTAR A DENSIDADE DE ELÉTRONS COM A ALTURA.

        Parameters
        ----------
        ne : PANDA SERIES - FLOATS
            DESCRIPTION.
        h : PANDA SERIES
            DESCRIPTION.
        data : STRING, optional
            DATA PARA O QUAL O DADO FOI ADQUIRIDO. The default is "".

        Returns
        -------
        None.

        """
        plt.figure(figsize=(5,5))
        plt.plot(ne,h,label = "$N_e$ "+data)
        plt.xlabel("Densidade de elétrons ($m^{-3}$)")
        plt.ylabel("Altura (km)")
       
        plt.title("Densidade de elétrons \n " + data)
        
        #plt.xlim(left=10e5)
        plt.legend()
        plt.grid()
        #$plt.xscale('log',base=10)
       
        plt.show()
    
        
    def plot_rhoion(self):
        plt.figure(figsize=(5,5))

        plt.plot(self.rhodado,self.dado["H(km)"])
        
        plt.ylabel("Altura (km)")
        plt.xlabel("Densidade de íons ($m^{-3}$)")
        plt.title("Densidade de ions \n ")
        
        plt.legend()
        plt.grid()
        plt.xscale("log")
        
        plt.show()


#nomearqIRI = "IRI2016_lat69_lon19_z300_s5_data30MAR2012_h6UT.txt"
#a = iri(nomearqIRI)
# #print(a.dado.dtypes)

# a.to_percentage()
# a.eletron_cm3to_m3()

# print(a.calc_rho_numion(a.dado["Ne(m3)"],a.rhodado["O+"]))
#print(a.dado["Te/K"])



# a.plot_densidade_e(a.dado["Ne(cm-3)"],a.dado["H(km)"],"21 de novembro de 2012")

# b = iri("IRI_lat-25_lon125_z0_750_s_5data01JAN2008_h12LT.txt")
# b.plot_densidade_e(b.dado["Ne(cm-3)"],b.dado["H(km)"]," 01 de JAN de 2008 às 12:00LT")

# c = iri("IRI_lat-25_lon125_z0_750_s_5data01JAN2008_h00LT.txt")
# c.plot_densidade_e(c.dado["Ne(cm-3)"],c.dado["H(km)"]," 01 de JAN de 2008 às 00:00LT")

# d = iri("IRI_lat-25_lon125_z0_750_s_5data01JAN2000_h00LT.txt")
# d.plot_densidade_e(d.dado["Ne(cm-3)"],d.dado["H(km)"]," 01 de JAN de 2000 às 00:00LT")

# e = iri("IRI_lat-25_lon125_z0_750_s_5data01JAN2000_h12LT.txt")
# e.plot_densidade_e(e.dado["Ne(cm-3)"],e.dado["H(km)"]," 01 de JAN de 2000 às 12:00LT")

# f = d.dado["Ne(cm-3)"]-c.dado["Ne(cm-3)"]
# d.plot_densidade_e(f,e.dado["H(km)"]," 2008 - 2000 às 12:00LT")

# g = b.dado["Ne(cm-3)"] - e.dado["Ne(cm-3)"]
# d.plot_densidade_e(g,b.dado["H(km)"]," 2008 - 2000 às 00:00LT")