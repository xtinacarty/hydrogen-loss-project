##LOAD IN DATA##
import xarray as xr; 

hydrogendata= input("name of Hydrogen Mixing Ratios data file: " )
h2odata = input("name of H2O Mixing Ratios data file: " )
tempatmodata = input("name of Atmospheric Temperature data file: " )
title = input("Choose a name for this run: ")
     


h_xr = xr.open_dataset(hydrogendata) #hydrogen data
ta_xr = xr.open_dataset(tempatmodata) #atmospheric temperature data
h2o_xr = xr.open_dataset(h2odata) #h2o data

#Pulling Values
    #Model Top is at lev[1] = 9.8e-6
ta_top = ta_xr.isel(lev=[1]) 
h_top= h_xr.isel(lev = [1]) 
h2o_top = h2o_xr.isel(lev = [1])
ta_therm = ta_xr.isel(lev = [3]) #Temp at ~ 100km altitude; lev[3] =2.67e-5

##SLICING OUT DAYSIDE VALUES##
    #day side = (90,270) deg longitude, substellar point at 180 deg lon
h_dayxr=h_top.sel(lon=slice(90,270))
h2o_dayxr=h2o_top.sel(lon=slice(90,270))
t_dayxr =ta_top.sel(lon=slice(90,270))
ttherm_dayxr = ta_therm.sel(lon=slice(90,270))

##SLICING OUT NIGHTSIDE VALUES##
    #night side = (0,90) and (270,360) deg longitude
h_nightxr = xr.concat([h_top.sel(lon=slice(0,90)), h_top.sel(lon=slice(270,360))], 'lon')
h2o_nightxr = xr.concat([h2o_top.sel(lon=slice(0,90)), h2o_top.sel(lon=slice(270,360))], 'lon')
t_nightxr = xr.concat([ta_top.sel(lon=slice(0,90)), ta_top.sel(lon=slice(270,360))], 'lon')
ttherm_nightxr = xr.concat([ta_therm.sel(lon=slice(0,90)), ta_therm.sel(lon=slice(270,360))], 'lon')

##DEFINING MEANS##
    #GLOBAL
Q = h_top.mean().H.values.tolist()
Qapprox = (h2o_top.mean().H2O.values.tolist())*2
T = ta_top.mean().T.values.tolist()
Ttherm = ta_therm.mean().T.values.tolist()

    #DAYSIDE
Q_day = h_dayxr.mean().H.values.tolist()
Qapprox_day = (h2o_dayxr.mean().H2O.values.tolist())*2
T_day = t_dayxr.mean().T.values.tolist()
Ttherm_day = ttherm_dayxr.mean().T.values.tolist()

    #NIGHTSIDE
Q_night = h_nightxr.mean().H.values.tolist()
Qapprox_night = (h2o_nightxr.mean().H2O.values.tolist())*2
T_night = t_nightxr.mean().T.values.tolist()
Ttherm_night = ttherm_nightxr.mean().T.values.tolist()
   


hatoms_total = (2.107 * 10**28)  #H atoms in the ocean (earth values)

###CALCULATIONS###
import hydrogenescape as fn

#Atmospheric Scale Height##
H = fn.atmoscaleheight(T) #global
H_day = fn.atmoscaleheight(T_day) #dayside
H_night = fn.atmoscaleheight(T_night) #nightside

##Escape Rates##
    #Global
r = fn.escaperate(Q, Ttherm, H)
rapprox = fn.escaperate(Qapprox, Ttherm, H)

    #Dayside
r_day = fn.escaperate(Q_day, Ttherm_day, H_day) 
rapprox_day = fn.escaperate(Qapprox_day, Ttherm_day, H_day) 

    #Nightside
r_night = fn.escaperate(Q_night, Ttherm_night, H_night)
rapprox_night = fn.escaperate(Qapprox_night, Ttherm_night, H_night)

##Ocean Survival Timescale##
hatoms = (2 * 10**29)  #H atoms cm**-2 in the ocean (earth values)

    #Global
gyrs = fn.oceansurvival(hatoms, r)
gyrs_approx = fn.oceansurvival(hatoms, rapprox)
    
    #Dayside
gyrs_day = fn.oceansurvival(hatoms*0.5, r_day)
gyrsapprox_day = fn.oceansurvival(hatoms*0.5, rapprox_day)
    #Nightside
gyrs_night = fn.oceansurvival(hatoms*0.5, r_night)
gyrsapprox_night = fn.oceansurvival(hatoms*0.5, rapprox_night)


##TABLE##
import plotly.graph_objects as go
from decimal import Decimal

        
values = [['<b>GLOBAL<b>', '<b>DAY SIDE<b>', '<b>NIGHT SIDE<b>'], #1st col
        ['%.2E' % Decimal(r), '%.2E' % Decimal(r_day), '%.2E' % Decimal(r_night)],
        ['%.2E' % Decimal(rapprox), '%.2E' % Decimal(rapprox_day), '%.2E' % Decimal(rapprox_night)],
        ['%.2E' % Decimal(gyrs), '%.2E' % Decimal(gyrs_day), '%.2E' % Decimal(gyrs_night)],
        ['%.2E' % Decimal(gyrs_approx), '%.2E' % Decimal(gyrsapprox_day), '%.2E' % Decimal(gyrsapprox_night)]]

fig = go.Figure(data=[go.Table(
  columnwidth = [125,250, 250, 250, 250],
  header = dict(
    values = [['RUN: %s'% title],
                  ['<b>H ESCAPE RATE</b> <br><i>(10<sup>7</sup> cm<sup>-2</sup> sec<sup>-1</sup>)</i>'],
             ['<b>H ESCAPE RATE<br>(APPROX.)</b><br><i>(10<sup>7</sup> cm<sup>-2</sup> sec<sup>-1</sup>)</i>'],
             ['<b>OCEAN SURVIVAL</b><br><i>(Billion Years)</i>'],
             ['<b>OCEAN SURVIVAL <br>(APPROX.) </b><br><i>(Billion Years)</i>']],
    line_color='white',
    fill_color=['#3B3561'],
    align=['left','left', 'left', 'left', 'left'],
    font=dict(color=['white'], size=[18, 14]),
    font_family=["courier new", "Abadi"],
    height=60,
  ),
  cells=dict(
    values=values,
    line_color='white',
    fill=dict(color=[['mistyrose', 'lemonchiffon', 'aliceblue']]),
    font_size=12,
      font_family = "Abadi",
    height=30)
    )
])
fig.show()

import plotly.io as pio
pio.kaleido.scope.default_width = 900

x = input("Export table?: [y/n] ")
if x == 'y':
    fig.write_image("Table(%s).png" % title)
    print('Table exported to "Table(%s).png"' % title)
elif x == 'n':
    pass
else:
    print("Invalid input-- table not exported")
    print("Input only accepts 'y' or 'n'")

##DICTIONARIES##
##DICTIONARIES##
if 'Variables' in locals():
    pass
else:
    Variables = {
        "Q": "Hydrogen Mixing Ratio",
        "Qapprox": "Hydrogen Mixing Ratio, approximated by H2O values",
        "T": "Model top temperature",
        "Ttherm": "Temperature of the Thermosphere, at approximately 100km altitude",
        "H": "Atmospheric scale height",
        "r": "Hydrogen Escape Rate",
        "rapprox": "Hydrogen Escape Rate using Qapprox",
        "gyrs": "Survival of 1 Earth Ocean per billion years given Q,",
        "gyrsapprox": "Survival of 1 Earth Ocean per billion years given Qapprox"
}
    
thisrun = (("Q_%s" % title, Q,), 
                   ("Qapprox_%s" % title, Qapprox), 
                   ("T_%s" % title, T),
                   ("Ttherm_%s" % title, Ttherm),
                   ("H_%s" % title, H), 
                   ("r_%s" % title, r), 
                   ("rapprox_%s" % title, rapprox), 
                   ("gyrs_%s" % title, gyrs), 
                   ("gyrs_approx_%s" % title, gyrs_approx))
if 'GlobalValues' in locals():
        GlobalValues.update(thisrun)
else:
    GlobalValues = {thisrun}
    
thisrunday = {
        "Q_%s" % title: Q_day,
        "Qapprox_%s" % title: Qapprox_day,
        "T_%s" % title: T_day,
        "Ttherm_%s" % title: Ttherm_day,
        "H_%s" % title: H,
        "r_%s" % title: r_day,
        "rapprox_%s" % title: rapprox_day,
        "gyrs_%s" % title: gyrs_day,
        "gyrsapprox_%s" % title: gyrsapprox_day}
if 'DaysideValues' in locals():
    DaysideValues.update(thisrunday)
else:
    DaysideValues = thisrunday

thisrunnight = {
        "Q_%s" % title: Q_night,
        "Qapprox_%s" % title: Qapprox_night,
        "T_%s" % title: T_night,
        "Ttherm_%s" % title: Ttherm_night,
        "H_%s" % title: H,
        "r_%s" % title: r_night,
        "rapprox_%s" % title: rapprox_night,
        "gyrs_%s" % title: gyrs_night,
        "gyrsapprox_%s" % title: gyrsapprox_night}
if 'NightsideeValues' in locals():
    NightsideValues.update(thisrunnight)
else:
    NightsideValues = thisrunnight


print("GlobalValues")
print(GlobalValues)
print("DaysideValues")
print(DaysideValues)
print("NightsideValues")
print(NightsideValues)
print("Variables")
print(Variables)
print("To run the program again, user can execute the following commands;")
print("$ import importlib")
print("$ importlib.reload(FinalProgram)")
    