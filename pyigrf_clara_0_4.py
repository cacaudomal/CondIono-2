#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:31:55 2022
Modificaçao de pyigrf
@author: clara
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pyIGRF: code to synthesise magnetic field values from the 13th generation of the
        International Geomagnetic Reference Field (IGRF), released in December 2019

 @author: Ciaran Beggan (British Geological Survey)
 
 See https://www.ngdc.noaa.gov/IAGA/vmod/ for information on the IGRF
 
 Based on existing codes: igrf13.f (FORTRAN) and chaosmagpy (Python3)
 
 With acknowledgements to: Clemens Kloss (DTU Space), David Kerridge (BGS),
      william Brown and Grace Cox.
 
     This is a program for synthesising geomagnetic field values from the 
     International Geomagnetic Reference Field series of models as agreed
     in December 2019 by IAGA Working Group V-MOD. 
     
     This is the 13th generation IGRF.
     
     The main-field models for 1900.0, 1905.0,..1940.0 and 2020.0 are 
     non-definitive, those for 1945.0, 1950.0,...2015.0 are definitive and
     the secular-variation model for 2020.0 to 2025.0 is non-definitive.

     Main-field models are to degree and order 10 (i.e. 120 coefficients)
     for 1900.0-1995.0 and to 13 (i.e. 195 coefficients) for 2000.0 onwards. 
     The predictive secular-variation model is to degree and order 8 (i.e. 80
     coefficients).

 Inputs:
 -------
     Inputs are via the command line:    
     
     Options include: 
         Write to (1) screen or (2) filename          
    
    Type of computation:
         (1) values at a single locations at one time (spot value)
         (2) values at same location at one year intervals (time series), 
         (3) grid of values at one time (grid); 
     
        Positions can be in:  
         (1) geodetic (WGS-84)
         (2) geocentric coordinates
         
     Latitude & longitude entered as: 
         (1) decimal degrees or 
         (2) degrees & minutes (not in grid)
    
     Altitude can be entered as:
         (1) Geodetic: height above the ellipsoid in km
         (2) Geocentric: distance from the Earth's centre in km
        
     Date: in decimal years (e.g. 2020.25)
 
 Outputs: 
 -----------
      : main field (in nanoTesla) and secular variation (in nanoTesla/year)
 
 
 Dependencies: 
 -------------
     : numpy, scipy
 
 Recent history of code:
 -----------------------
     Initial release: April 2020 (Ciaran Beggan, BGS)
 
"""

from scipy import interpolate
import igrf_utils as iut
import io_options_clara as ioo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
#import pyglow 
# import os.path
# from os import path

# Load in the file of coefficients
IGRF_FILE = r'./IGRF13.shc'
igrf = iut.load_shcfile(IGRF_FILE, None)

def string_para_float(dado):
    """
    Método para separar um string contendo só números em uma lista de floats.
    
    Parameters
    ----------
    dado : VETOR DE STRING
        DADOS A SEREM CONVERTIDOS.

    Returns
    -------
    dado : VETOR DE FLOAT
        Dado convertido para float.  
        
    """
    for x in range(len(dado)):
    	dado[x] = float(dado[x])
    #print("\nDado:\n",dado)
    return dado

class IGRF():
    """
    CLASS with functions to facilitate the usage of the  Ciaran Beggan's 
    2019 library:
    pyIGRF: code to synthesise magnetic field values from the 13th generation of the
            International Geomagnetic Reference Field (IGRF), released in December 2019.
    
    INFO:
    Df grid is a DATAFRAME of the shape:
        ['Latitude','Longitude','Altitude','Year','Declination' , 'Inclination' , 'Horizontal_intensity' , 'Total_intensity','North_component','East_component','Vertical_component' , 'DeclinationSV' , 'InclinationSV' , 'HorizontalSV' , 'TotalSV' , 'NorthSV' , 'EastSV' , 'VerticalSV',"B(T)"]
    in which :
        
    Declination (D) em graus
    Inclination (I) em graus
    Horizontal intensity (H) nT
    Total intensity (F)      nT
    North component (X)      nT
    East component (Y)       nT
    Vertical component (Z)   nT
    Declination SV (D)       arcmin/yr
    Inclination SV (I)       arcmin/yr
    Horizontal SV (H)        nT/yr
    Total SV (F)             nT/yr
    North SV (X)             nT/yr
    East SV (Y)              nT/yr
    Vertical SV (Z)          nT/yr
    
    ...
    
     Dependencies: 
     -------------
         : numpy, scipy,pyigrf, matplotlib, pandas, geopandas
         
        
    Attributes
    ----------
        lat : FLOAT
            Latitude inicial em graus decimais.
        lon : FLOAT
            LONGITUDE INICIAL EM GRAUS DECIMAIS.
        h : FLOAT
            ALTURA INICIAL [km].
        ano : INT
            ANO NO QUAL SERÁ CALCULADO OS VALORES DO CAMPO GEOMAGNÉTICO.
        name_saida : STRING
            NOME DO ARQUIVO DE SAÍDA SEM O FORMATO.
    
    Methods
    ----------
    
        set_name_saida(self, name_saida):
        _salva_dataframe(self,nomesaida,ResultStore)
        get_grid(self,nomearq)
        _calc_igrf(self,parLat,parLon,parH,parAno)
        calc_grid(self,intervalo_h = 10,lim_h = 500,intervalo_lat=-5,lim_lat=-60,intervalo_lon=-5,lim_lon=-110)
        calc_perfilh(self,lim_h = 300,intervalo_h = 5)
        coord_centro(self,dado,nomearq = "centro_coord")
        create_matriz_IntensidadeTotal(self,DfGrid,h = 400)
        calc_gradiente(self,DfGrid)
        plot_grid(self,Dfgrid,h = 400)
        _nT_to_T_grid(self)
        _nT_to_T_perfil
        
    """
    
    def __init__(self,lat,lon,h,ano,name_saida):
        """
        Inicializa a classe IGRF.
        
        Parameters
        ----------
        lat : FLOAT
            Latitude inicial em graus decimais.
        lon : FLOAT
            LONGITUDE INICIAL EM GRAUS DECIMAIS.
        h : FLOAT
            ALTURA INICIAL [km].
        ano : INT
            ANO NO QUAL SERÁ CALCULADO OS VALORES DO CAMPO GEOMAGNÉTICO.
        name_saida : STRING
            NOME DO ARQUIVO DE SAÍDA SEM O FORMATO.

        Returns
        -------
        None.

        """
        self.set_name_saida(name_saida)        
        self.entrada_usuario = (lat,lon,h,ano)
    
    def set_name_saida(self, name_saida):
        """
        SETS THE OUTPUT FILE NAME. 

        Parameters
        ----------
        name_saida : STRING
            NAME OF THE OUTPUT PILE TO BE USED.

        Returns
        -------
        None.

        """
        self.name_saida = name_saida
    
    def _salva_dataframe(self,nomesaida,ResultStore):
        """
        SALVA OS DADOS GUARDADOS NUM DATAFRAME NUM ARQUIVO.
        
        parameters
        -------
        nomesaida : STRING
            NOME DO ARQUIVO NO QUAL SERÃO SALVOS OS DADOS SEM O FORMATO.
            
        ResultStore : DATAFRAME
            VALORES CALCULADOS QUE SERÃO SALVOS NO ARQUIVO. 

        """
        
        ResultStore.to_csv(nomesaida+".csv",sep=",",header=True)
        print("pyigrf_clara_0_3 - _salva_dataframe : arquivo", nomesaida," aberto.\n")
    
    def get_grid(self,nomearq):
        """
        CRIA UM DATAFRAME A PARTIR DE UM ARQUIVO CSV COM VALORES PREVIAMENTE CALCULADOS.

        Parameters
        ----------
        nomearq : STRING
            NOME DO ARQUIVO ONDE ESTÃO OS DADOS.

        Returns
        -------
        self.Dfgrid : DATAFRAME
            DATAFRAME COM O GRID COM TODOS OS VALORES CALCULADOS.

        """        
        
        self.Dfgrid = pd.read_csv(nomearq,sep=",")
        return self.Dfgrid
    
    def _calc_igrf(self,parLat,parLon,parH,parAno):
        """
        CALCULA AS COMPONENTES DO CAMPO MAGNÉTICO EM UM ÚNICO PONTO.
        Calcula no ponto escolhido os seguinte componentes do campo geomagnético: 
        Declination (D) em graus
        Inclination (I) em graus
        Horizontal intensity (H) nT
        Total intensity (F)      nT
        North component (X)      nT
        East component (Y)       nT
        Vertical component (Z)   nT
        Declination SV (D)       arcmin/yr
        Inclination SV (I)       arcmin/yr
        Horizontal SV (H)        nT/yr
        Total SV (F)             nT/yr
        North SV (X)             nT/yr
        East SV (Y)              nT/yr
        Vertical SV (Z)          nT/yr
        
        Parameters
        ----------
        name : STRING
            nome do arquivo de entrada
        dado : LIST
            Lista com os parametros para configurar os dados que vc quer do programa
            [opção, lat, long, altura, ano]
        parLat : 
            
        parLon : 
            
        parH : 
            
        parAno :
        
        Returns
        -------
        eff: FLOAT
            Total intensity  in [nT]
        date : FLOAT
            ano decimal
        alt : FLOAT
            altitude [km]
        lat : FLOAT
            Latitude in degree
        lon : FLOAT
            longitude in degree
        dec : FLOAT
            Declination in degree  
        inc : FLOAT
            Inclination in degree
        hoz : FLOAT
            Horizontal intensity [nT]
        eff : FLOAT
            Total intensity [nT]
        X : FLOAT
            North component [nT]
        Y : FLOAT
            East component  [nT]
        Z : FLOAT
            Vertical component [nT]
        decs : FLOAT
            Declination SV   arcmin/yr
        incs : FLOAT
            Inclination SV   arcmin/yr
        hozs : FLOAT
            Horizontal SV  [nT/yr]
        effs : FLOAT
            Total SV (F) [nT/yr]
        dX : FLOAT
            North SV [nT/yr]
        dY : FLOAT
            East SV  [nT/yr]
        dZ : FLOAT
            Vertical SV  [nT/yr]
       
        """        
        #print("\n\n IGRF nfunc_calc_igrf parametros: ", self.parametros)
        #print( '\nCalculating values at one location and date')
        date, alt, lat, colat, lon, itype, sd, cd = ioo.option1(1,parLat,parLon,parH,parAno)
    
            
        # Interpolate the geomagnetic coefficients to the desired date(s)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        f = interpolate.interp1d(igrf.time, igrf.coeffs, fill_value='extrapolate')
        coeffs = f(date)    
        
        # Compute the main field B_r, B_theta and B_phi value for the location(s) 
        Br, Bt, Bp = iut.synth_values(coeffs.T, alt, colat, lon,
                                  igrf.parameters['nmax'])
        
        # For the SV, find the 5 year period in which the date lies and compute
        # the SV within that period. IGRF has constant SV between each 5 year period
        # We don't need to subtract 1900 but it makes it clearer:
        epoch = (date-1900)//5    
        epoch_start = epoch * 5
        
        # Add 1900 back on plus 1 year to account for SV in nT per year (nT/yr):
        coeffs_sv = f(1900+epoch_start+1) - f(1900+epoch_start)   
        Brs, Bts, Bps = iut.synth_values(coeffs_sv.T, alt, colat, lon,
                                  igrf.parameters['nmax'])
        
        # Use the main field coefficients from the start of each five epoch
        # to compute the SV for Dec, Inc, Hor and Total Field (F) 
        # [Note: these are non-linear components of X, Y and Z so treat separately]
        coeffsm = f(1900+epoch_start);
        Brm, Btm, Bpm = iut.synth_values(coeffsm.T, alt, colat, lon,
                                  igrf.parameters['nmax'])
        
        # Rearrange to X, Y, Z components 
        X = -Bt; Y = Bp; Z = -Br
        # For the SV
        dX = -Bts; dY = Bps; dZ = -Brs 
        Xm = -Btm; Ym = Bpm; Zm = -Brm
        
        # Rotate back to geodetic coords if needed
        if (itype == 1):
            t = X; X = X*cd + Z*sd;  Z = Z*cd - t*sd
            t = dX; dX = dX*cd + dZ*sd;  dZ = dZ*cd - t*sd
            t = Xm; Xm = Xm*cd + Zm*sd;  Zm = Zm*cd - t*sd
            
        # Compute the four non-linear components 
        dec, hoz, inc, eff = iut.xyz2dhif(X,Y,Z)
        
        # The IGRF SV coefficients are relative to the main field components 
        # at the start of each five year epoch e.g. 2010, 2015, 2020
        decs, hozs, incs, effs = iut.xyz2dhif_sv(Xm, Ym, Zm, dX, dY, dZ)
        
        #print("\n\nfunc_main: eff: ",eff,"\n")
            
        #Colocando nas coordenadas e altura certos
        alt, lat = iut.geo_to_gg(alt, colat)
        lat = 90-lat
        
        return lat,lon, np.round(alt, decimals=3), date, dec, inc, hoz, eff, X, Y, Z, decs, hozs,\
            incs, effs, dX, dY, dZ
                     
    def calc_grid(self,intervalo_h = 10,lim_h = 500,intervalo_lat=-5,lim_lat=-60,intervalo_lon=-5,lim_lon=-110):
        """
        CRIA O GRID PARA OS INTERVALOS DE ALTURA, LATITUDE E LONGITUDE ESPECIFICADOS E SALVA OS DADOS NUM ARQUIVO. 
        
        Parameters
        ----------
        intervalo_h : FLOAT, optional
            espaçamento da altura entre os pontos calculado [km]. The default is 5.
        lim_h : FLOAT, optional
            ALTURA MÁXIMA [km]. The default is 300.
        intervalo_lat : FLOAT, optional
            INTERVALO ENTRE OS PONTOS DE LATITUDE A SEREM CALCULADOS [decimal degree]. The default is -1.
        lim_lat : FLOAT, optional
            LATITUDE FINAL [decimal degree]. The default is -35.
        intervalo_lon : FLOAT, optional
            INTERVALO ENTRE OS PONTOS DE LONGITUDE A SEREM CALCULADOS [decimal degree]. The default is -1.
        lim_lon : FLOAT, optional
            LONGITUDE MÁXIMA [decimal degree]. The default is -35.

        Returns
        -------
        lgrid : DATAFRAME
            DATAFRAME COM OS VALORES CALCULADOS NAS ALTURAS E COORDENADAS DADAS.

        """
        lgrid = []
        print("pyigrf_0_3 - calc_grid self.Dfgrid : starting to calculate grid")
        for lat in np.arange(self.entrada_usuario[0],lim_lat,intervalo_lat):
            #print("lat",lat)
            for lon in np.arange(self.entrada_usuario[1],lim_lon,intervalo_lon):
                #print("lon",lon)
                for h in np.arange(self.entrada_usuario[2],lim_h,intervalo_h): 
                    lgrid.append(self._calc_igrf(lat,lon,h,self.entrada_usuario[3])) #quando ele calcula ele já guarda no arquivo o valor daquele ponto
                    #print("\npyigrf_clara_0_2  calc_igrf lat :",lat)
                    #print("calc_grid",lat,lon,h,lgrid)
                    
        #print("\n\ncalc_grid - grid :",len(grid))
        self.Dfgrid = pd.DataFrame(lgrid,columns = ['Latitude','Longitude','Altitude','Year','Declination' , 'Inclination' , 'Horizontal_intensity' , 'Total_intensity','North_component','East_component','Vertical_component' , 'DeclinationSV' , 'InclinationSV' , 'HorizontalSV' , 'TotalSV' , 'NorthSV' , 'EastSV' , 'VerticalSV'] )
        self._nT_to_T_grid()
        #print("pyigrf_0_3 - calc_grid self.Dfgrid head : ", self.Dfgrid.head())
        
        nomearq = self.name_saida + "_grid" 
        self._salva_dataframe(nomearq,self.Dfgrid)
        
        print("pyigrf_clara - calc_grid : o grid foi salvo")
        
        return lgrid
        
    def calc_perfilh(self,lat,lon,lim_h = 750,intervalo_h = 5):
        """
        Função que calcula os valores do campo magnético em um perfil de altitude.

        Parameters
        ----------
        lim_h : INT, optional
            LIMITE SUPERIOR PARA O QUAL A ALTURA SERÁ CALCULADA. The default is 300.
        intervalo_h : INT, optional
            INTERVALOS PARA O QUAL OS VALORES SERÃO CALCULADOS. The default is 5.

        Returns
        -------
        self.DfPerfilh : DATAFRAME
            PERFIL DE ALTURA PARA OS VALORES DO CAMPO MAGNÉTICO EM LATITUDE E LONGITUDE FIXA.
        """
        nomesaida = self.name_saida + "_perfil"
        lperfil = []
        print("pyigrf_clara - calc_perfilh: Calculando Perfil de Altura do Campo Magnético")
        for i in np.arange(self.entrada_usuario[2],lim_h,intervalo_h): 
            lperfil.append(self._calc_igrf(lat,lon,i,self.entrada_usuario[3])) #quando ele calcula ele já guarda no arquivo o valor daquele ponto
        
        #print("pyigrf_clara - calc_perfilh : lperfil",lperfil[1])
        
        self.DfPerfilh = pd.DataFrame(lperfil,columns =['Latitude','Longitude','Altitude','Year','Declination' , 'Inclination' , 'Horizontal_intensity' , 'Total_intensity','North_component','East_component','Vertical_component' , 'DeclinationSV' , 'InclinationSV' , 'HorizontalSV' , 'TotalSV' , 'NorthSV' , 'EastSV' , 'VerticalSV'] )
        self._nT_to_T_perfil()
        
        self._salva_dataframe(nomesaida,self.DfPerfilh) # Salva os valores num arquivo
        
        return self.DfPerfilh
    
    def coord_centro(self,dado,nomearq = "centro_coord"):
        """
        Identifica onde está o centro da anomalia magnética do atlântico sul e em que altura. 
        Considera o centro como o ponto de menor valor de intensidade do campo independentemente da altura.

        Parameters
        ----------
        dado : DATAFRAME
            DADOS CALCULADOS PELA FUNÇÃO CALC_GRID.

        Returns
        -------
       linha_dataframe["Latitude"] : PANDA SERIES
            DATA SERIES COM OS VALORES DE LATITUDE EM GRAUS DECIMAIS.
        
        linha_dataframe["Longitude"] : PANDA SERIES
            DATA SERIES COM OS VALORES DE LONGITUDE EM GRAUS DECIMAIS.
        
        linha_dataframe["Total_intensity"] : PANDA SERIES
            DATA SERIES COM VALOR DA INTENSIDADE DO CAMPO MAGNÉTICO EM [nT]

        """
        #valorB_centro = dado["Total_intensity"].min()
        #print("coord_centro Tipo dado:",type(dado))
        
        posicao_linha_dataframe = dado["Total_intensity"].idxmin() # qual linha esta o menor valor? 
        print("coord_centro posicao linha dataframe : " ,posicao_linha_dataframe,"len dado",len(dado))
        
        linha_dataframe = dado.iloc[posicao_linha_dataframe] #pega essa linha e salva no data frame
        self._salva_dataframe(self.name_saida + "_coord_centro", linha_dataframe)
        
        return linha_dataframe["Latitude"],linha_dataframe["Longitude"], linha_dataframe["Total_intensity"],linha_dataframe["Altitude"]
    
    def coord_centro_all_heights(self,dado,nomearq = "centro_coord"):
        """
        Identifica onde está o centro da anomalia magnética do atlântico sul. 
        Considera o centro como o ponto de menor valor de intensidade do campo.

        Parameters
        ----------
        dado : DATAFRAME
            DADOS CALCULADOS PELA FUNÇÃO CALC_GRID.

        Returns
        -------
       linha_dataframe["Latitude"] : DATA SERIES
            DATA SERIES COM OS VALORES DE LATITUDE EM GRAUS DECIMAIS.
        
        linha_dataframe["Longitude"] : DATA SERIES
            DATA SERIES COM OS VALORES DE LONGITUDE EM GRAUS DECIMAIS.
        
        linha_dataframe["Total_intensity"] : DATA SERIES
            DATA SERIES COM VALOR DA INTENSIDADE DO CAMPO MAGNÉTICO EM [nT]

        """
        #valorB_centro = dado["Total_intensity"].min()
        #print("coord_centro Tipo dado:",type(dado))
        
        "*** EM PROGRESSO ***" 
        
        posicao_linha_dataframe = dado["Total_intensity"].idxmin() # qual linha esta o menor valor? 
        #print("coord_centro posicao linha dataframe : " ,posicao_linha_dataframe,"len dado",len(dado))
        
        linha_dataframe = dado.iloc[posicao_linha_dataframe] #pega essa linha e salva no data frame
        self._salva_dataframe(self.name_saida + "_coord_centro", linha_dataframe)
        
        return linha_dataframe["Latitude"],linha_dataframe["Longitude"], linha_dataframe["Total_intensity"],linha_dataframe["Altitude"]

    def plot_grid(self,Dfgrid,h = 400):
        """
        PLOTA A INTENSIDADE DO CAMPO MAGNÉTICO NUMA DADA ALTURA SOBRE O MAPA DA TERRA.

        Parameters
        ----------
        Dfgrid : DATAFRAME
            DATAFRAME COM O GRID DOS VALORES CALCULADOS NA FUNÇÃO calc_grid, TEM QUE TER UMA COLUNA "Total_intensity".
            TEM QUE SER UM GRID QUADRADO NA LAT E LONG. 
        h : INTEGER, optional
            ALTURA PARA A QUAL O CAMPO MAGNÉTICO SERÁ PLOTADO [km]. The default is 400.

        Returns
        -------
        None.

        """
        a = Dfgrid.loc[Dfgrid['Altitude'] == h] #dados só para altitude h
        
        # PLOTA OS PAISES NO MAPA
        fig, ax = plt.subplots(figsize=(8,8))
        countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
        countries.plot(ax = ax, color = "lightgrey")
        countries.head()
        
        #embelezando o mapa
        ax.set_title("Mapa da Intensidade de B ano de 2020")
        ax.set_ylabel("Latitude")
        ax.set_xlabel("Longitude")
        #ax.legend()
        ax.grid(visible=True, alpha=0.5)
        
        #f, ax = plt.subplots(1,2, sharex=True, sharey=True)
        #plota mapa de cores
        ax.tripcolor(a["Longitude"],a["Latitude"],a["Total_intensity"]) 
        
        #plota o mapa de contornos
        #nota: mudar cor das linhas
        ax.tricontour(a["Longitude"],a["Latitude"],a["Total_intensity"],15)
        
        plt.show()
    
    def _nT_to_T_grid(self):
        self.Dfgrid["B(T)"] = self.Dfgrid["Total_intensity"] * 10**-9
    
    def _nT_to_T_perfil(self):
        self.DfPerfilh["B(T)"] = self.DfPerfilh["Total_intensity"] * 10**-9
        
##======================================================
#
#PROGRAMA PRINCIPAL
#
#========================================
#name_saida = "laat" 

#dado = IGRF(0,50,80,2008,name_saida)
#dado.calc_grid(intervalo_h = 10,lim_h = 300,intervalo_lat=-10,lim_lat=-60,intervalo_lon=-10,lim_lon=-120)
#hmax = 500 #a altura maxima para o qual o igrf vai ser calculado
# #print("dado criado")

#dado.calc_grid(lim_h = 500)
#dado.plot_grid(dado.Dfgrid,400)

#h = 400 #altura para a qual vou calcular o centro da anomalia

#c = dado.Dfgrid["Altitude"] - dado.Dfgrid["Latitude"]
#print ("\n\nC: ",c)

#grid_h = dado.Dfgrid.loc[dado.Dfgrid['Altitude'] == h]


#print("Dados da coordenada do ponto central",dado.coord_centro(grid_h,dado.name_saida))

#dado.create_matriz_IntensidadeTotal(dado.Dfgrid)
#f = dado.calc_gradiente(dado.Dfgrid)

#plt.tripcolor(dado.Dfgrid["Longitude"],dado.Dfgrid["Latitude"],dado.Dfgrid["Total_intensity"]) #plota mapa de cores

#plotando o campo de vetores do gradiente
#axis 0 é a latitude  (que vai ser o y)
#axis 1 é a longitude (que é o x)
#plt.quiver(dado.Dfgrid["Longitude"].unique(),dado.Dfgrid["Latitude"].round(6).unique(),f[1],f[0], linewidth=None, color='pink')

#a = dado.calc_grid(intervalo_h = 10,lim_h=hmax,intervalo_lat=-10,lim_lat=-60,intervalo_lon=-10,lim_lon=-90)

#calculando a norma dos vetores para plotar
#criando o vetor para receber os valores calculados
#o eixo x é latitude e o y é longitude
#norma = np.zeros((len(dado.Dfgrid["Latitude"].round(6).unique()),len(dado.Dfgrid["Longitude"].unique())))

#===== Testes Perfil =====================
#dado.calc_perfilh(hmax,10)
# print(dado.DfPerfilh)
# print("\n\n",dado.DfPerfilh.columns.values)
#plt.plot(dado.DfPerfilh['Total_intensity'],dado.DfPerfilh['Altitude'])

#===== Plot Grid =============================

# a = dado.grid.loc[dado.grid['Altitude'] == 400] #dados só para altitude 400
# b = dado.dado.loc[dado.dado['Altitude']==250] #dados só para altitude 250
# #c=a.loc[['Lat']==0]

# #print(len(c['Lat']))
# #a.plot(x="Lon", y="Lat", kind="scatter", c="Total_intensity", colormap="YlOrRd",title="Scatter plot total intensity B" )

# print(len(dado.grid["Latitude"].round(6).unique()),len(dado.grid["Longitude"].unique()))
# g = f[0] + f[1]
# #plt.tricontour(dado.grid["Latitude"].round(6).unique(),dado.grid["Longitude"].unique(),g,15)
# plt.plot(g[1],dado.grid["Latitude"].round(6).unique())
# plt.plot(g[0],dado.grid["Longitude"].round(6).unique())
# #print(val.reshape(5,5))

# #contourf([a['Lon'],a['Lat'],] a["Total_intensity"], , **kwargs)
# #plt.imshow(a[['Lat','Lon','Total_intensity']],extent=[-30,-35, -30, -35],cmap ="RdYlBu",interpolation='bilinear')
# plt.show()   