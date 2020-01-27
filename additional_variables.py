#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 21 19:53:02 2019

@author: research
"""


import pandas as pd
import numpy as np

pathto_raw_data = '...'

#----------------------------------------------------------------------------------------------------------------------------------

pathto_1_200 = '...'
pathto_200_400 = '...'
pathto_400_500 = '...'
pathto_500_600 = '...'

#----------------------------------------------------------------------------------------------------------------------------------

pathtoaddvars_1_200 = '...'
pathtoaddvars_200_400 = '...'

#----------------------------------------------------------------------------------------------------------------------------------

quarterlycompvars = pd.read_csv( ( pathto_raw_data + 'quartcomp.csv' ), usecols=['LPERMNO','cusip','GVKEY','datadate','cik','sic','naics','ibq','saleq','txtq','revtq','cogsq','prccq','cshoq','xsgaq','atq','actq','cheq','lctq','dlcq','ppentq','ceqq','seqq','pstkq','atq','ltq','pstkrq', 'dlttq', 'xrdq', 'txpq', 'dpq', 'oancfy', 'invtq', 'spiq', 'rectq', 'xintq', 'ppegtq', 'capxy', 'gdwlq', 'apq', 'obkq', 'acoq', 'intanq', 'aoq', 'lcoq', 'loq', 'drcq', 'niq', 'scstkcy'])
quarterlycompvars.drop_duplicates(subset=['datadate', 'LPERMNO'], inplace=True)
quarterlycompvars.set_index('datadate', inplace=True)
quarterlycompvars.index = pd.to_datetime(quarterlycompvars.index, format='%d/%m/%Y').to_period('M').asfreq('D')

#----------------------------------------------------------------------------------------------------------------------------------
''' Creating Daily Variables'''

# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'RET']).to_csv('./data/WRDSdata/092019/13092019/dailyret.csv')
#
# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'PRC']).to_csv('./data/WRDSdata/092019/13092019/dailyprc.csv')
#
# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'VOL']).to_csv('./data/WRDSdata/092019/13092019/dailyvol.csv')
#
# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'SHROUT']).to_csv('./data/WRDSdata/092019/13092019/dailyshrout.csv')
#
# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'ASKHI']).to_csv('./data/WRDSdata/092019/13092019/dailyaskhi.csv')
#
# pd.read_csv(daily_crsp, usecols=['PERMNO','date', 'BIDLO']).to_csv('./data/WRDSdata/092019/13092019/dailybidlo.csv')

#----------------------------------------------------------------------------------------------------------------------------------

'''CRSP MONTHLY INFORMATION'''

monthlycrspvars = pd.read_csv((pathto_raw_data + 'monthlycrsp.csv'), usecols=['PERMNO','date', 'SHROUT', 'ALTPRC', 'CUSIP', 'EXCHCD', 'SHRCD', 'SICCD', 'SPREAD','VOL'])

monthlycrspvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
monthlycrspvars.set_index('date', inplace=True)
monthlycrspvars.index = pd.to_datetime(monthlycrspvars.index, format='%d/%m/%Y').to_period('M').asfreq('D')

#----------------------------------------------------------------------------------------------------------------------------------

    # Market Cap

        # Calculating Monthly Market Cap

monthlycrspvars['mve_f'] = monthlycrspvars['ALTPRC'].abs() * monthlycrspvars['SHROUT']
monthlycrspvars['mve_f'].replace(0, np.nan, inplace=True)

#----------------------------------------------------------------------------------------------------------------------------------

    # Compress cusip 1:6
quarterlycompvars['cnum'] = quarterlycompvars['cusip'].apply(lambda x: str(x)[:8])

    # Substr sic 1:2
quarterlycompvars['sic2'] = quarterlycompvars['sic'].apply(lambda x: str(x)[:2])
monthlycrspvars['SICCD2'] = monthlycrspvars['SICCD'].apply(lambda x: str(x)[:2])

#----------------------------------------------------------------------------------------------------------------------------------
    # Count
quarterlycompvars['count'] = 1
monthlycrspvars['count'] = 1

quarterlycompvars['count'] = pd.DataFrame(quarterlycompvars.groupby(by='GVKEY').cumsum()['count'])
monthlycrspvars['count'] = pd.DataFrame(monthlycrspvars.groupby(by='PERMNO').cumsum()['count'])

#----------------------------------------------------------------------------------------------------------------------------------

'''						*do some clean up, several of these variables have lots of missing values;

                        if missing(xint) then xint0=0;
                            else xint0=xint;
                        if missing(xsga) then xsga0=0;
                            else xsga0=0;
                        
                            '''

    # xint0 calculation
quarterlycompvars['xintq0'] = quarterlycompvars['xintq'].replace(np.nan,0)

    # xsga0 calculation
quarterlycompvars['xsgaq0'] = quarterlycompvars['xsgaq'].replace(np.nan,0)

'''
# =============================================================================
# DOWNGRADED ANOMALIES
# =============================================================================
'''
#----------------------------------------------------------------------------------------------------------------------------------
 
    # 9 bm ( Book to Market )
'''book-to-Market = (Common Equity / Market Cap)
'''
bmvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
bmvars = bmvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

bm = bmvars / mve_f

bm.astype(float).to_csv(pathto_200_400 +'209.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 33 ep ( Earnings to Price )
'''Earning-to-Price = ( Income Before Extra Ordinary Items / Market Cap )
'''
epvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
epvars = epvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

ep = epvars / mve_f

ep.astype(float).to_csv(pathto_200_400 +'233.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 13 cashpr ( Cash Productivity )
'''Cash Productivity = ( ( Market Cap + Long Term Debt - Total Asset) / Cash and Short Term Investments)
'''
ltdebtvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlttq')
ltdebtvars = ltdebtvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)


mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

cashpr = ((mve_f + ltdebtvars - atvars) / chevars.replace(0,np.nan))

cashpr.astype(float).to_csv(pathto_200_400 +'213.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''TOTAL SUM OF DIVIDENDS IS NOT BEING FOUND IN THE QUARTERLY DATABASE
#    # 30 dy (Dividend to Price)
#dy = dvt / mve_f
#
#dvtvars = quarterlycompvars.pivot(columns='LPERMNO', values='dvt')
#dvtvars = dvtvars.fillna(value=None, method='ffill', limit=3)
#
#mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')
#
#dy = dvtvars / mve_f
#
#dy.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/30.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 43 lev ( Leverage )
'''lev = lt / mve_f
'''
levvars = quarterlycompvars.pivot(columns='LPERMNO', values='ltq')
levvars = levvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

lev = levvars / mve_f

lev.astype(float).to_csv(pathto_200_400 +'243.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 86 sp ( Sales to Price )
'''sp = sale / mve_f
'''
spvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
spvars = spvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

sp = spvars / mve_f

sp.astype(float).to_csv(pathto_200_400 +'286.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' EBIT IS NOT BEING FOUND IN THE QUARTERLY DATABASED
#    # 77 roic ( Return on Ivested Capital )
#roic = (ebit - nopi) / (ceq + lt - che)
#
#ebitvars = quarterlycompvars.pivot(columns='LPERMNO', values='ebit')
#ebitvars = ebitvars.fillna(value=None, method='ffill', limit=3)
#nopivars = quarterlycompvars.pivot(columns='LPERMNO', values='nopi')
#nopivars = nopivars.fillna(value=None, method='ffill', limit=3)
#ceqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
#ceqvars = ceqvars.fillna(value=None, method='ffill', limit=3)
#ltvars = quarterlycompvars.pivot(columns='LPERMNO', values='ltq')
#ltvars = ltvars.fillna(value=None, method='ffill', limit=3)
#chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
#chevars = chevars.fillna(value=None, method='ffill', limit=3)
#
#denominator = ceqvars + ltvars - chevars
#denominator = denominator.replace(0,np.nan)
#
#roic = (ebitvars - nopivars) / denominator
#
#roic.astype(float).to_csv(pathto_200_400 +'277.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 71 rd_sale ( R&D to sale )
'''rd_sale = xrd / sale
'''
xrdvars = quarterlycompvars.pivot(columns='LPERMNO', values='xrdq')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=3)
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)

rd_sale = xrdvars / salevars.replace(0,np.nan)

rd_sale.astype(float).to_csv(pathto_200_400 +'271.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 70 rd_mve ( R&D to Market Capitalization )
'''rd_mve = xrd / mve_f
'''
xrdvars = quarterlycompvars.pivot(columns='LPERMNO', values='xrdq')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

rd_mve = xrdvars / mve_f

rd_mve.astype(float).to_csv(pathto_200_400 +'270.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 5 agr ( Asset Growth )
'''agr = (at / lag(at)) - 1
'''
agrvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
agrvars = agrvars.fillna(value=None, method='ffill', limit=3)

agr = (agrvars / agrvars.replace(0,np.nan).shift(3)) - 1

agr.astype(float).to_csv(pathto_200_400 +'205.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 34 gma ( Gross Profitability )
'''gma = (revt - cogs) / lag(at)
'''
revtvars = quarterlycompvars.pivot(columns='LPERMNO', values='revtq')
revtvars = revtvars.fillna(value=None, method='ffill', limit=3)
cogsvars = quarterlycompvars.pivot(columns='LPERMNO', values='cogsq')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

gma = (revtvars - cogsvars) / atvars.replace(0,np.nan).shift(3)

gma.astype(float).to_csv(pathto_200_400 +'234.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''THIS IS ALREADY MONTHLY DATA

#    # 17 chcsho ( Change in Shares Outstanding )
#chcsho = (csho / lag(csho)) - 1
#
#csho = monthlycrspvars.pivot(columns='PERMNO', values='SHROUT')
#
#chcsho = (csho / csho.replace(0,np.nan).shift()) - 1
#
#chcsho.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/17.csv')
#
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 44 lgr ( Growth in Long Term Debt )
'''lgr = (lt / lag(lt)) - 1'''
lgrvars = quarterlycompvars.pivot(columns='LPERMNO', values='ltq')
lgrvars = lgrvars.fillna(value=None, method='ffill', limit=3)

lgr = (lgrvars / lgrvars.replace(0,np.nan).shift(3)) - 1

lgr.astype(float).to_csv(pathto_200_400 +'244.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 2 acc ( Accruals )
'''	acc=(ib-oancf) /  ((at+lag(at))/2);
if missing(oancf) then acc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/  ((at+lag(at))/2)
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

        # Creating mask
accmask = (oancfvars.isna())

        # Calculating First Variable
acc = (ibvars - oancfvars) / ((atvars + atvars.shift(3)).replace(0,np.nan) / 2)

        # Applying mask
actvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actvars = actvars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
lctvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctvars = lctvars.fillna(value=None, method='ffill', limit=3)
dlcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlcq')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=3)
txpvars = quarterlycompvars.pivot(columns='LPERMNO', values='txpq')
txpvars = txpvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

acc[accmask] =( ( ( (actvars - actvars.shift(3)) -
                             (chevars - chevars.shift(3)) ) -
                           ( (lctvars - lctvars.shift(3)) -
                             (dlcvars - dlcvars.shift(3)) -
                             (txpvars - txpvars.shift(3)) -
                             (dpvars) ) ) /
                           ( (atvars + atvars.shift(3)).replace(0,np.nan) / 2) )
    
acc.astype(float).to_csv(pathto_200_400 +'202.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 65 pctacc ( Percent Accruals )
'''pctacc=(ib-oancf)/abs(ib);
if ib=0 then pctacc=(ib-oancf)/.01;
if missing(oancf) then pctacc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/abs(ib);
if missing(oancf) and ib=0 then pctacc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/.01;
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)

        # Creating Masks
pctaccmask1 = (ibvars==0)
pctaccmask2 = (oancfvars.isna())
pctaccmask3 = ((oancfvars.isna()) & (ibvars==0))

        # Calculating first variable
pctacc = (ibvars - oancfvars) / ibvars.abs()

        # Applying First Mask
pctacc[pctaccmask1] = (ibvars - oancfvars) / .01

        # Applying Second Mask
actvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actvars = actvars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
lctvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctvars = lctvars.fillna(value=None, method='ffill', limit=3)
dlcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlcq')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=3)
txpvars = quarterlycompvars.pivot(columns='LPERMNO', values='txpq')
txpvars = txpvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)

pctacc[pctaccmask2] = ( ( ( (actvars - actvars.shift(3)) -
                            (chevars - chevars.shift(3)) ) -
                          ( (lctvars - lctvars.shift(3)) -
                            (dlcvars - dlcvars.shift(3)) -
                            (txpvars - txpvars.shift(3)) -
                            (dpvars) ) ) /
                          ( (ibvars.abs()) ) )

        # Applying Third Mask
pctacc[pctaccmask3] = ( ( ( (actvars - actvars.shift(3)) -
                            (chevars - chevars.shift(3)) ) -
                          ( (lctvars - lctvars.shift(3)) -
                            (dlcvars - dlcvars.shift(3)) -
                            (txpvars - txpvars.shift(3)) -
                            (dpvars) ) ) /
                          ( (.01) ) )

pctacc.astype(float).to_csv(pathto_200_400 +'265.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # 14 cfp ( Cash Flow to Price Ratio )
'''cfp=oancf/mve_f
if missing(oancf) then cfp=(ib-(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) ))/mve_f
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)
actvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actvars = actvars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
lctvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctvars = lctvars.fillna(value=None, method='ffill', limit=3)
dlcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlcq')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=3)
txpvars = quarterlycompvars.pivot(columns='LPERMNO', values='txpq')
txpvars = txpvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

        # Creating mask
cfpmask = (oancfvars.isna())

        # Calculating First Variable
cfp = oancfvars / mve_f

        # Applying mask
cfp[cfpmask] = ( (ibvars -
    ( ( (actvars - actvars.shift(3)) -
        (chevars - chevars.shift(3)) ) -
      ( (lctvars - lctvars.shift(3)) -
        (dlcvars - dlcvars.shift(3)) -
        (txpvars - txpvars.shift(3)) -
        (dpvars) ) ) ) /
         mve_f )


cfp.astype(float).to_csv(pathto_200_400 + '214.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # 1 absacc ( Absolute Accruals )
'''Absolute of Accruals
Use Variable #2 acc ( Accruals)
'''
acc = pd.read_csv(pathto_200_400 + '202.csv')
acc.set_index('datadate', inplace=True)
acc.index = pd.to_datetime(acc.index).to_period('M').asfreq('D')

absacc = acc.abs()

absacc.astype(float).to_csv(pathto_200_400 + '201.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 4 age ( # of years Since First Compustat Coverage )
'''	if first.gvkey then count=1;
	else count+1;
Original Paper's Description: Firm Age is defined as the number of months between event month t and the first month that a stock appears in CRSP'''
agevars = quarterlycompvars.pivot(columns='LPERMNO', values='count')
agevars = agevars.fillna(value=None, method='ffill', limit=3)

age = agevars.copy()

age.astype(float).to_csv(pathto_200_400 + '204.csv')

#----------------------------------------------------------------------------------------------------------------------------------
	# 19 chinv ( Change in Inventory )
'''chinv=(invt-lag(invt))/((at+lag(at))/2)
'''
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

chinv = ( invtvars - invtvars.shift(3) ) / ( ( atvars + atvars.shift(3)).replace(0,np.nan) / 2 )

chinv.astype(float).to_csv(pathto_200_400 + '219.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 1 spii
'''if spi ne 0 and not missing(spi) then spii=1; else spii=0
'''
quarterlycompvars['spiiq'] = 0
quarterlycompvars['spiiq'][(quarterlycompvars['spiq']!=0) & (quarterlycompvars['spiq'].notnull())] = 1

spii = quarterlycompvars.pivot(columns='LPERMNO', values='spiiq')
spii = spii.fillna(value=None, method='ffill', limit=3)

spii.astype(float).to_csv(pathtoaddvars_200_400 + '201.csv')


#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 2 spi ( Altered spi )
'''spi=spi/((at+lag(at))/2)
'''
spivars = quarterlycompvars.pivot(columns='LPERMNO', values='spiiq')
spivars = spivars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

spi = spivars / ( (atvars + atvars.shift(3)) / 2 ).replace(0,np.nan)

spi.astype(float).to_csv(pathtoaddvars_200_400 + '202.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 3 cf
'''cf=oancf/((at+lag(at))/2);
if missing(oancf) then cf=(ib-(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) ))/((at+lag(at))/2)
'''
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

        # Creating mask
cfmask = (oancfvars.isna())

        # Calculating First Variable
cf = oancfvars / ( (atvars + atvars.shift(3)) / 2 ).replace(0,np.nan)

        # Applying Mask
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
actvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actvars = actvars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
lctvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctvars = lctvars.fillna(value=None, method='ffill', limit=3)
dlcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlcq')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=3)
txpvars = quarterlycompvars.pivot(columns='LPERMNO', values='txpq')
txpvars = txpvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)

cf[cfmask] = ( (ibvars -
    ( ( (actvars - actvars.shift(3)) -
        (chevars - chevars.shift(3)) ) -
      ( (lctvars - lctvars.shift(3)) -
        (dlcvars - dlcvars.shift(3)) -
        (txpvars - txpvars.shift(3)) -
        (dpvars) ) ) ) /
      ( (atvars + atvars.shift(3)) / 2 ).replace(0,np.nan) )

cf.astype(float).to_csv(pathtoaddvars_200_400 + '203.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''TOTAL EMPLOYEES IS NOT BEING FOUND IN THE QUARTERLY DATABASE
#	# 38 hire ( Employee Growth Rate )
#hire=(emp-lag(emp))/lag(emp);
#if missing(emp) or missing(lag(emp)) then hire=0
#
#empvars = quarterlycompvars.pivot(columns='LPERMNO', values='emp')
#empvars = empvars.fillna(value=None, method='ffill', limit=3)
#
#hire = (empvars - empvars.shift(3)) / empvars.shift(3)
#hire[(empvars.isna()) & (empvars.shift(3)).isna()] = 0
#hire = hire.replace(np.inf, 0)
#
#hire = hire.replace(-np.inf, 0)
#
#hire.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/38.csv')
'''

#----------------------------------------------------------------------------------------------------------------------------------
	# 84 ( Sales Growth )
'''sgr=(sale/lag(sale))-1
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)

sgr = ( salevars / salevars.shift(3).replace(0,np.nan) ) - 1

sgr.astype(float).to_csv(pathto_200_400 + '284.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 4 chpm
'''chpm = (ib / sale) - (lag(ib) / lag(sale))
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)

chpm = ( ibvars / salevars.replace(0,np.nan) ) - ( ibvars.shift(3) / salevars.shift(3).replace(0,np.nan) )

chpm.astype(float).to_csv(pathtoaddvars_200_400 + '204.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 5 chato
'''chato=(sale/((at+lag(at))/2)) - (lag(sale)/((lag(at)+lag2(at))/2));
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

chato=( salevars / ((atvars + atvars.shift(3)) / 2).replace(0,np.nan) ) - ( salevars.shift(3) / ((atvars.shift(3) + atvars.shift(6)) / 2).replace(0,np.nan) )

chato.astype(float).to_csv(pathtoaddvars_200_400 + '205.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 61 pchsale_pchinvt ( % Change in Sales - % Change in Inventory )
'''pchsale_pchinvt=((sale-lag(sale))/lag(sale))-((invt-lag(invt))/lag(invt))
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

pchsale_pchinvt=( (salevars - salevars.shift(3) ) / salevars.shift(3).replace(0,np.nan) ) - ( (invtvars - invtvars.shift(3) ) / invtvars.shift(3).replace(0,np.nan) )

pchsale_pchinvt.astype(float).to_csv(pathto_200_400 + '261.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 62 pchsale_pchrect ( % Change in Sales - % Change in A/R )
'''pchsale_pchrect=((sale-lag(sale))/lag(sale))-((rect-lag(rect))/lag(rect))
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
rectvars = quarterlycompvars.pivot(columns='LPERMNO', values='rectq')
rectvars = rectvars.fillna(value=None, method='ffill', limit=3)

pchsale_pchrect = ( (salevars - salevars.shift(3) ) / salevars.shift(3).replace(0,np.nan) ) - ( (rectvars - rectvars.shift(3) ) / rectvars.shift(3).replace(0,np.nan) )

pchsale_pchrect.astype(float).to_csv(pathto_200_400 + '262.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 59 pchgm_pchsale ( % Change in Gross Margin - % Change in Sales )
'''pchgm_pchsale=(((sale-cogs)-(lag(sale)-lag(cogs)))/(lag(sale)-lag(cogs)))-((sale-lag(sale))/lag(sale))
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
cogsvars = quarterlycompvars.pivot(columns='LPERMNO', values='cogsq')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=3)

pchgm_pchsale = ( ( ( salevars - cogsvars ) - ( salevars.shift(3) - cogsvars.shift(3) ) ) / ( salevars.shift(3) - cogsvars.shift(3) ).replace(0,np.nan) ) - ( ( salevars - salevars.shift(3) ) / salevars.shift(3).replace(0,np.nan) )

pchgm_pchsale.astype(float).to_csv(pathto_200_400 + '259.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 63 pchsale_pchxsga ( % Change in Sales - % Change in SG&A )
'''pchsale_pchxsga=( (sale-lag(sale))/lag(sale) )-( (xsga-lag(xsga)) /lag(xsga) )
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
xsgavars = quarterlycompvars.pivot(columns='LPERMNO', values='xsgaq0')
xsgavars = xsgavars.fillna(value=None, method='ffill', limit=3)

pchsale_pchxsga = ( (salevars - salevars.shift(3) ) / salevars.shift(3).replace(0,np.nan) ) - ( (xsgavars - xsgavars.shift(3) ) / xsgavars.shift(3).replace(0,np.nan) )

pchsale_pchxsga.astype(float).to_csv(pathto_200_400 + '263.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 26 depr ( Depreciation / PP&E )
'''depr=dp/ppent
'''
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)

depr = dpvars / ppentvars.replace(0,np.nan)

depr.astype(float).to_csv(pathto_200_400 + '226.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 58 pchdepr ( % Change in Depreciation )
'''pchdepr=((dp/ppent)-(lag(dp)/lag(ppent)))/(lag(dp)/lag(ppent))
'''
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)

pchdepr = ( ( dpvars / ppentvars.replace(0,np.nan) ) - ( dpvars.shift(3) / ppentvars.shift(3).replace(0,np.nan) ) ) / ( dpvars.shift(3) / ppentvars.shift(3).replace(0,np.nan) ).replace(0,np.nan)

pchdepr.astype(float).to_csv(pathto_200_400 + '258.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''XAD IS NOT BEING FOUND IN THE QUARTERLY DATABASE
#    # Additional Variable 6 chadv
#chadv=log(1+xad)-log((1+lag(xad)));	*had error here before, might work better now...
#
#xadvars = quarterlycompvars.pivot(columns='LPERMNO', values='xad')
#xadvars = xadvars.fillna(value=None, method='ffill', limit=3)
#
#chadv = np.log( xadvars + 1 ) - np.log( xadvars.shift(3) + 1 )
#
#chadv.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/6.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 42 invest ( Capital Expenditures and Inventory )
'''invest=( 	(ppegt-lag(ppegt)) +  (invt-lag(invt))	)	/ lag(at)
if missing(ppegt) then invest=( 	(ppent-lag(ppent)) +  (invt-lag(invt))	)	/ lag(at)
'''
ppegtvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppegtq')
ppegtvars = ppegtvars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
investmask = (ppegtvars.isna())

        # Calculating First Variable
invest = ( ( ppegtvars - ppegtvars.shift(3) ) + ( invtvars - invtvars.shift(3) ) ) / atvars.shift(3).replace(0,np.nan)

        # Applying Mask
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)

invest[investmask] = ( ( ppentvars - ppentvars.shift(3) ) + ( invtvars - invtvars.shift(3) ) ) / atvars.shift(3).replace(0,np.nan)

invest.astype(float).to_csv(pathto_200_400 + '242.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 32 egr ( Growth in Common Shareholder Equity
'''egr=( (ceq-lag(ceq))/lag(ceq)  )
'''
ceqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=3)

egr = ( ceqvars - ceqvars.shift(3) ) / ceqvars.shift(3).replace(0,np.nan)

egr.astype(float).to_csv(pathto_200_400 + '232.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variables 7 capx
'''if missing(capx) and count>=2 then capx=ppent-lag(ppent)
Use count variable
'''
capxvars = quarterlycompvars.pivot(columns='LPERMNO', values='capxy')
capxvars = capxvars.fillna(value=None, method='ffill', limit=3)
countvars = quarterlycompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=3)
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
capxmask = ( (capxvars.isna()) & (countvars >= 2) )

        # Applying Mask
capxvars[capxmask] = ppentvars - ppentvars.shift(3)

capx = capxvars.copy()

capx.astype(float).to_csv(pathtoaddvars_200_400 + '207.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 8 pchcapx
'''pchcapx = (capx - lag(capx)) / lag(capx)
Use additional Variable 7 ( capx )
'''
pchcapx = ( capx - capx.shift(3) ) / capx.shift(3).replace(0,np.nan)

pchcapx.astype(float).to_csv(pathtoaddvars_200_400 + '208.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 35 grCAPX ( Growth in Capital Expenditures )
'''grcapx = (capx - lag2(capx)) / lag2(capx)
Use additional Variable 7 ( capx )
'''
grcapx = ( capx - capx.shift(6) ) / capx.shift(6).replace(0,np.nan)

grcapx.astype(float).to_csv(pathto_200_400 + '235.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 9 grGW
'''grGW = (gdwl - lag(gdwl)) / lag(gdwl);
if missing(gdwl) or gdwl=0 then grGW=0;
if gdwl ne 0 and not missing(gdwl) and missing(grGW) then grGW=1;
'''
quarterlycompvars['grGWmask'] = np.nan

gdwlvars = quarterlycompvars.pivot(columns='LPERMNO', values='gdwlq')
gdwlvars = gdwlvars.fillna(value=None, method='ffill', limit=3)

        # Calculating First Variable
grGW = ( gdwlvars - gdwlvars.shift(3 )) / gdwlvars.shift(3).replace(0,np.nan)

        # Creating Mask
grgwmask1 = ( (gdwlvars.isna()) | (gdwlvars==0) )
grgwmask2 = ( (gdwlvars!=0) & (gdwlvars.notnull()) & (grGW.isna()) )

        # Applying Mask
grGW[grgwmask1] = 0
grGW[grgwmask2] = 1

grGW.astype(float).to_csv(pathtoaddvars_200_400 + '209.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''GWO IS NOT BEING FOUND IN THE QUARTERLY DATABASE
#    # Additional Variable 10 woGW
#if (not missing(gdwlia) and gdwlia ne 0) or (not missing(gdwlip) and gdwlip ne 0) or (not missing(gwo) and gwo ne 0) then woGW=1;
#else woGW=0;
#
#gdwliavars = quarterlycompvars.pivot(columns='LPERMNO', values='gdwlia')
#gdwliavars = gdwliavars.fillna(value=None, method='ffill', limit=3)
#gdwlipvars = quarterlycompvars.pivot(columns='LPERMNO', values='gdwlip')
#gdwlipvars = gdwlipvars.fillna(value=None, method='ffill', limit=3)
#gwovars = quarterlycompvars.pivot(columns='LPERMNO', values='gwo')
#gwovars = gwovars.fillna(value=None, method='ffill', limit=3)
#
#        # Calculating First Variable
#woGW = pd.DataFrame(columns=gdwliavars.columns, index=gdwliavars.index, data=0)
#
#        # Creating Mask
#wogwmask = ( ( (gdwliavars.notnull()) & (gdwliavars!=0) ) | ( (gdwlipvars.notnull()) & (gdwlipvars!=0) ) | ( (gwovars.notnull()) & (gwovars!=0) ) )
#
#        # Applying Mask
#woGW[wogwmask] = 1
#
#woGW.astype(float).to_csv(pathtoaddvars_200_400 + '210.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 91 tang ( Dept Capacity / Firm Tangibility )
'''tang = (che + rect * 0.715 + invt * 0.547 + ppent * 0.535) / at
'''
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
rectvars = quarterlycompvars.pivot(columns='LPERMNO', values='rectq')
rectvars = rectvars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

tang = (chevars + rectvars * 0.715 + invtvars * 0.547 + ppentvars * 0.535) / atvars.replace(0,np.nan)

tang.astype(float).to_csv(pathto_200_400 + '291.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 85 sin ( Sin Stocks )
'''if (2100 <= sic <= 2199) or (2080 <= sic <= 2085) or (
naics in ('7132', '71312', '713210', '71329', '713290', '72112', '721120'))
then
sin = 1; else sin = 0;
'''
        # Creating Mask
sinmask = ( ( ( 2100 <= quarterlycompvars['sic']) & (quarterlycompvars['sic'] <= 2199) ) | ( (2080 <= quarterlycompvars['sic']) & (quarterlycompvars['sic'] <= 2085) ) | ( quarterlycompvars['naics'].isin(['7132', '71312', '713210', '71329', '713290', '72112', '721120']) ) )

        # Calculating First Variable
quarterlycompvars['sin'] = 0

        # Applying Mask
quarterlycompvars['sin'][sinmask] = 1

sin = quarterlycompvars.pivot(columns='LPERMNO', values='sin')
sin = sin.fillna(value=None, method='ffill', limit=3)

sin.astype(float).to_csv(pathto_200_400 + '285.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 11 act
'''if missing(act) then act=che+rect+invt
'''
actvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actvars = actvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
actmask = ( actvars.isna() )

        # Calculating First Variable
act = actvars.copy()

        # Applying Mask
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)
rectvars = quarterlycompvars.pivot(columns='LPERMNO', values='rectq')
rectvars = rectvars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

act[actmask] = chevars + rectvars + invtvars

act.astype(float).to_csv(pathtoaddvars_200_400 + '211.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 12 lct
'''if missing(lct) then lct=ap
'''
lctvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctvars = lctvars.fillna(value=None, method='ffill', limit=3)
apvars = quarterlycompvars.pivot(columns='LPERMNO', values='apq')
apvars = apvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
lctmask = ( lctvars.isna() )

        # Calculating First Variable
lct = lctvars.copy()

        # Applying Mask
lct[lctmask] = apvars

lct.astype(float).to_csv(pathtoaddvars_200_400 + '212.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 25 currat ( Current Ratio )
'''currat = act / lct
Use Additional Variable 11 ( act ) 879
Use Additional Variable 12 ( lct ) 902
'''

currat = act / lct.replace(0,np.nan)

currat.astype(float).to_csv(pathto_200_400 + '225.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 57 pchcurrat ( % Change in Current Ratio )
'''pchcurrat = ((act / lct) - (lag(act) / lag(lct))) / (lag(act) / lag(lct))
Use Additional Variable 11 ( act ) 879
Use Additional Variable 12 ( lct ) 902
'''
pchcurrat = ( ( act / lct.replace(0,np.nan) ) - ( act.shift(3) / lct.shift(3).replace(0,np.nan) ) ) / ( act.shift(3) / lct.shift(3).replace(0,np.nan) ).replace(0,np.nan)

pchcurrat.astype(float).to_csv(pathto_200_400 + '257.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 68 quick ( Quick Ratio )
'''quick = (act - invt) / lct
Use Additional Variable 11 ( act ) 879 
Use Additional Variable 12 ( lct ) 902
'''
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

quick = ( act - invtvars.astype(float) ) / lct.replace(0,np.nan)

quick.astype(float).to_csv(pathto_200_400 + '268.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 60 pchquick ( % Change in quick Ratio )
'''pchquick = (	(act -invt ) /lct - (lag(act ) -lag(invt) ) /lag(lct) )/  (   (lag(act) - lag(invt)) / lag(lct))
Use Additional Variable 9 ( act ) 879
Use Additional Variable 10 ( lct ) 902
'''
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

pchquick = (	(act - invtvars ) / lct.replace(0,np.nan) - ( act.shift(3)  - invtvars.shift(3) ) / lct.shift(3).replace(0,np.nan) ) /  (  ( act.shift(3) -  invtvars.shift(3) ) / lct.shift(3) ).replace(0,np.nan)

pchquick.astype(float).to_csv(pathto_200_400 + '260.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 79 salecash ( Sales to Cash )
'''salecash = sale / che
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
chevars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
chevars = chevars.fillna(value=None, method='ffill', limit=3)

salecash = salevars / chevars.replace(0,np.nan)

salecash.astype(float).to_csv(pathto_200_400 + '279.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 81 salerec ( Sales to Receivables )
'''salerec = sale / rect
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
rectvars = quarterlycompvars.pivot(columns='LPERMNO', values='rectq')
rectvars = rectvars.fillna(value=None, method='ffill', limit=3)

salerec = salevars / rectvars.replace(0,np.nan)

salerec.astype(float).to_csv(pathto_200_400 + '281.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 80 saleinv ( Sales to Inventory )
'''saleinv = sale / invt
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

saleinv = salevars / invtvars.replace(0,np.nan)

saleinv.astype(float).to_csv(pathto_200_400 + '280.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 64 pchsaleinv ( % Change Sales-to-Inventory )
'''pchsaleinv = ((sale / invt) - (lag(sale) / lag(invt))) / (lag(sale) / lag(invt))
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)

pchsaleinv = ( ( salevars / invtvars.replace(0,np.nan) ) - ( salevars.shift(3) / invtvars.shift(3).replace(0,np.nan) ) ) / ( salevars.shift(3) / invtvars.shift(3).replace(0,np.nan) ).replace(0,np.nan)

pchsaleinv.astype(float).to_csv(pathto_200_400 + '264.csv')


#----------------------------------------------------------------------------------------------------------------------------------
    # 12 cashdebt ( Cash Flow to Debt )
'''cashdebt = (ib + dp) / ((lt + lag(lt)) / 2)
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)
ltvars = quarterlycompvars.pivot(columns='LPERMNO', values='ltq')
ltvars = ltvars.fillna(value=None, method='ffill', limit=3)

cashdebt =  ( ibvars + dpvars ) / ( ( ltvars + ltvars.shift(3) )  / 2 ).replace(0,np.nan)

cashdebt.astype(float).to_csv(pathto_200_400 + '212.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' FATB, FATL IS NOT BEING FOUND IN THE QUARTERLY DATABASE
# 
#    # 72 realestate ( Real Estate Holdings )
#realestate = (fatb + fatl) / ppegt
#if missing(ppegt) then realestate=(fatb+fatl) / ppent
#
#fatbvars = quarterlycompvars.pivot(columns='LPERMNO', values='fatbq')
#fatbvars = fatbvars.fillna(value=None, method='ffill', limit=3)
#fatlvars = quarterlycompvars.pivot(columns='LPERMNO', values='fatlq')
#fatlvars = fatlvars.fillna(value=None, method='ffill', limit=3)
#ppegtvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppegtq')
#ppegtvars = ppegtvars.fillna(value=None, method='ffill', limit=3)
#ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
#ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)
#
#fatbvars[fatbvars[fatlvars==0].isna()] = 0
#
#        # Creating Mask
#realestatemask = (ppegtvars.isna())
#
#        # Calculating First Variable
#realestate = ( fatbvars + fatlvars) / ppegtvars
#
#        # Applying Mask
#realestate[realestatemask] = ( fatbvars + fatlvars) / ppentvars.replace(0,np.nan)
#
#realestate.astype(float).to_csv(pathto_200_400 + '272.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
''' TOTAL DIVIDENS IS NOT BEING FOUND IN THE QUARTERLY DATABASE

    # 27 divi ( Dividend Initiation )
if (not missing(dvt) and dvt > 0) and (lag(dvt)=0 or missing(lag(dvt))) then divi=1; else divi=0

dvtvars = quarterlycompvars.pivot(columns='LPERMNO', values='dvt')
dvtvars = dvtvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
divimask = ( ( (dvtvars.notnull()) & (dvtvars>0) ) & (dvtvars.shift(3)==0 ) | ( dvtvars.shift(3).isna() ) )

        # Calculating First Variable
divi = pd.DataFrame(columns=dvtvars.columns, index=dvtvars.index, data=0)

        # Applying Mask
divi[divimask] = 1

divi.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/27.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
''' TOTAL DIVIDENS IS NOT BEING FOUND IN THE QUARTERLY DATABASE

    # 28 divo ( Dividen Omission )
if (missing(dvt) or dvt=0) and (lag(dvt) > 0 and not missing(lag(dvt))) then divo=1; else divo=0

dvtvars = quarterlycompvars.pivot(columns='LPERMNO', values='dvt')
dvtvars = dvtvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
divomask = ( ( (dvtvars.isna()) | (dvtvars==0) ) & (dvtvars.shift(3)>0 ) & ( dvtvars.shift(3).notnull() ) )

        # Calculating First Variable
divo = pd.DataFrame(columns=dvtvars.columns, index=dvtvars.index, data=0)

        # Applying Mask
divo[divomask] = 1

divo.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/28.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 13 obklg
'''obklg = ob / ((at + lag(at)) / 2)
'''
obvars = quarterlycompvars.pivot(columns='LPERMNO', values='obkq')
obvars = obvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

obklg = obvars / ( atvars + atvars.shift(3) / 2 ).replace(0,np.nan)

obklg.astype(float).to_csv(pathtoaddvars_200_400 + '213.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 14 chobklg
'''chobklg = (ob - lag(ob)) / ((at + lag(at)) / 2)
'''
obvars = quarterlycompvars.pivot(columns='LPERMNO', values='obkq')
obvars = obvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

chobklg = ( obvars - obvars.shift(3) ) / ( atvars + atvars.shift(3) / 2 ).replace(0,np.nan)

chobklg.astype(float).to_csv(pathtoaddvars_200_400 + '214.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' DM IS NOT BEING FOUND IN THE QUARTERLY DATABASE
    # 83 secureind ( Secured Debt Indicator )
if not missing(dm) and dm ne 0 then securedind=1; else securedind=0

dmvars = quarterlycompvars.pivot(columns='LPERMNO', values='dm')
dmvars = dmvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
secureindmask = ( (dmvars.notnull()) & (dmvars!=0) )

        # Calculating First Variable
secureind = pd.DataFrame(columns=dmvars.columns, index=dmvars.index, data=0)

        # Applying Mask
secureind[secureindmask] = 1

secureind.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/83.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
''' DM, DLTT IS NOT BEING FOUND IN THE QUARTERLY DATABASE
    # 82 secured ( Secured Debt )
secured = dm / dltt

dmvars = quarterlycompvars.pivot(columns='LPERMNO', values='dm')
dmvars = dmvars.fillna(value=None, method='ffill', limit=3)
dlttvars = quarterlycompvars.pivot(columns='LPERMNO', values='dltt')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=3)

secured = dmvars / dlttvars.replace(0,np.nan)

secured.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/82.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
''' DC, CSHRC IS NOT BEING FOUND IN THE QUARTERLY DATABASE
    # 24 convind ( convertible Debt Indicator )
if not missing(dc) and dc ne 0 or ( not missing(cshrc) and CSHRC ne 0) then convind=1; else convind=0

dcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dc')
dcvars = dcvars.fillna(value=None, method='ffill', limit=3)
cshrcvars = quarterlycompvars.pivot(columns='LPERMNO', values='cshrc')
chhrcvars = cshrcvars.fillna(value=None, method='ffill', limit=3)

    # Creating Mask
convindmask = ( ( (dcvars.notnull()) & (dcvars!=0) ) | ( (cshrcvars.notnull()) & (cshrcvars!=0) ) )

    # Calculating First Variable
convind = pd.DataFrame(columns=dcvars.columns, index=dcvars.index, data=0)

    # Applying Mask
convind[convindmask] = 1

convind.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/24.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
''' DC IS NOT BEING FOUND IN THE QUARTERLY DATABASE
    # Additional Variable 15 conv
conv = dc / dltt

dcvars = quarterlycompvars.pivot(columns='LPERMNO', values='dc')
dcvars = dcvars.fillna(value=None, method='ffill', limit=3)
dlttvars = quarterlycompvars.pivot(columns='LPERMNO', values='dltt')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=3)

conv = dcvars / dlttvars.replace(0,np.nan)

conv.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/15.csv')
'''

#----------------------------------------------------------------------------------------------------------------------------------

    # 36 grltnoa ( Growth in Long Term Net Operating Assets )
'''grltnoa = ((rect + invt + ppent + aco + intan + ao - ap - lco - lo) - 
(lag(rect) + lag(invt) + lag(ppent) + lag(aco) + lag(intan) + lag(ao) - lag(ap) - lag(lco) - lag(lo))- 
(rect - lag(rect) + invt - lag(invt) + aco - lag(aco) - (ap - lag(ap) + lco - lag(lco)) - dp)) / 
((at + lag(at)) / 2)
'''
rectvars = quarterlycompvars.pivot(columns='LPERMNO', values='rectq')
rectvars = rectvars.fillna(value=None, method='ffill', limit=3)
invtvars = quarterlycompvars.pivot(columns='LPERMNO', values='invtq')
invtvars = invtvars.fillna(value=None, method='ffill', limit=3)
ppentvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=3)
acovars = quarterlycompvars.pivot(columns='LPERMNO', values='acoq')
acovars = acovars.fillna(value=None, method='ffill', limit=3)
intanvars = quarterlycompvars.pivot(columns='LPERMNO', values='intanq')
intanvars = intanvars.fillna(value=None, method='ffill', limit=3)
aovars = quarterlycompvars.pivot(columns='LPERMNO', values='aoq')
aovars = aovars.fillna(value=None, method='ffill', limit=3)
apvars = quarterlycompvars.pivot(columns='LPERMNO', values='apq')
apvars = apvars.fillna(value=None, method='ffill', limit=3)
lcovars = quarterlycompvars.pivot(columns='LPERMNO', values='lcoq')
lcovars = lcovars.fillna(value=None, method='ffill', limit=3)
lovars = quarterlycompvars.pivot(columns='LPERMNO', values='loq')
lovars = lovars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)

grltnoa = ( ( ( (rectvars + invtvars + ppentvars + acovars + intanvars + aovars - apvars - lcovars - lovars) -
( rectvars.shift(3) + invtvars.shift(3) + ppentvars.shift(3) + acovars.shift(3) + intanvars.shift(3) + aovars.shift(3) - apvars.shift(3) - lcovars.shift(3) - lovars.shift(3) ) ) -
( rectvars - rectvars.shift(3) + invtvars - invtvars.shift(3) + acovars - acovars.shift(3) - ( apvars - apvars.shift(3) + lcovars - lcovars.shift(3) ) - dpvars ) ) /
( (atvars + atvars.shift(3) ) / 2).replace(0,np.nan) )

grltnoa.astype(float).to_csv(pathto_200_400 + '/236.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 16 chdrc
'''chdrc = (dr - lag(dr)) / ((at + lag(at)) / 2)
'''
drvars = quarterlycompvars.pivot(columns='LPERMNO', values='drcq')
drvars = drvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

chdrc = ( drvars - drvars.shift(3) ) / ( ( atvars + atvars.shift(3) ) / 2).replace(0,np.nan)

chdrc.astype(float).to_csv(pathtoaddvars_200_400 + '216.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 69 rd ( R&D Increase )
'''if ((xrd / at) - (lag(xrd / lag(at)))) / (lag(xrd / lag(at))) > .05 then rd=1; else rd=0
'''
xrdvars = quarterlycompvars.pivot(columns='LPERMNO', values='xrdq')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)
countvars = quarterlycompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
rdmask = ( (( (xrdvars / atvars.replace(0,np.nan)) - ( ( xrdvars / atvars.shift(3).replace(0,np.nan) ) ) ) / ( xrdvars / atvars.shift(3).replace(0,np.nan) )) > 0.5 )

        # Calculating First Variable
rd = pd.DataFrame(columns=xrdvars.columns, index=xrdvars.index, data=0)

        # Applying Mask
rd[rdmask] = 1

countvars[countvars.notnull()] = 1
rd *=countvars

rd.astype(float).to_csv(pathto_200_400 + '269.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variables 17 rdbias
'''rdbias = (xrd / lag(xrd)) - 1 - ib / lag(ceq)
'''
xrdvars = quarterlycompvars.pivot(columns='LPERMNO', values='xrdq')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=3)
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
ceqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=3)

rdbias = (xrdvars / xrdvars.shift(3).replace(0,np.nan)) - 1 - (ibvars / ceqvars.shift(3).replace(0,np.nan))

rdbias.astype(float).to_csv(pathtoaddvars_200_400 + '217.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 18 roe
'''roe = ib / lag(ceq)
'''
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
ceqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=3)

roe = ibvars / ceqvars.shift(3).replace(0,np.nan)

roe.astype(float).to_csv(pathtoaddvars_200_400 + '218.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 54 operprof ( Operating Profitability )
'''operprof = (revt-cogs-xsga0-xint0)/lag(ceq)
'''
revtvars = quarterlycompvars.pivot(columns='LPERMNO', values='revtq')
revtvars = revtvars.fillna(value=None, method='ffill', limit=3)
cogsvars = quarterlycompvars.pivot(columns='LPERMNO', values='cogsq')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=3)
xsga0vars = quarterlycompvars.pivot(columns='LPERMNO', values='xsgaq0')
xsga0vars = xsga0vars.fillna(value=None, method='ffill', limit=3)
xint0vars = quarterlycompvars.pivot(columns='LPERMNO', values='xintq0')
xint0vars = xint0vars.fillna(value=None, method='ffill', limit=3)
ceqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=3)

operprof = (revtvars-cogsvars-xsga0vars-xint0vars) / ceqvars.shift(3).replace(0,np.nan)

operprof.astype(float).to_csv(pathto_200_400 + '254.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 67 ps ( Financial Statement Score )
'''ps = (ni>0)+(oancf>0)+(ni/at > lag(ni)/lag(at))+(oancf>ni)+
(dltt/at < lag(dltt)/lag(at))+(act/lct > lag(act)/lag(lct))+
((sale-cogs)/sale > (lag(sale)-lag(cogs))/lag(sale))+ 
(sale/at > lag(sale)/lag(at))+ (scstkc=0)
Use Additional Variable 11 ( act )
Use Additional Variable 12 ( lct )
'''
actvars = act.copy()
lctvars = lct.copy()
nivars = quarterlycompvars.pivot(columns='LPERMNO', values='niq')
nivars = nivars.fillna(value=None, method='ffill', limit=3)
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)
dlttvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlttq')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=3)
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)
cogsvars = quarterlycompvars.pivot(columns='LPERMNO', values='cogsq')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=3)
scstkcvars = quarterlycompvars.pivot(columns='LPERMNO', values='scstkcy')
scstkcvars = scstkcvars.fillna(value=None, method='ffill', limit=3)

ps = ( ( nivars [ nivars > 0 ] ) +
       ( oancfvars [ oancfvars > 0 ] ) +
       ( ( nivars / atvars.replace(0,np.nan) ) [ ( nivars / atvars.replace(0,np.nan) ) > ( nivars.shift(3) / atvars.shift(3).replace(0,np.nan) ) ] ) +
       ( oancfvars [ oancfvars > nivars ] ) +
       ( ( dlttvars / atvars.replace(0,np.nan) ) [ ( dlttvars / atvars.replace(0,np.nan) ) < ( dlttvars.shift(3) / atvars.shift(3).replace(0,np.nan) ) ] ) +
       ( ( actvars / lctvars.replace(0,np.nan) ) [ ( actvars / lctvars.replace(0,np.nan) ) > ( actvars.shift(3) / lctvars.shift(3).replace(0,np.nan) ) ] ) +
       ( ( ( salevars - cogsvars ) / salevars ) [ ( ( salevars - cogsvars ) / salevars ) > ( ( salevars.shift(3) - cogsvars.shift(3) ) / salevars.shift(3).replace(0,np.nan) )  ] ) +
       ( ( salevars / atvars.replace(0,np.nan) ) [ ( salevars / atvars )  > ( salevars.shift(3) / atvars.shift(3) ) ] ) +
       ( scstkcvars [ scstkcvars == 0 ] ) )

ps.astype(float).to_csv(pathto_200_400 + '267.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''TXFO and TXFED IS NOT BEING FOUND IN THE QUARTERLY DATABASE
# Additional Variable 19 tb_1
*-----Lev and Nissim(2004);
if fyear <= 1978 then tr=.48;
if 1979 <= fyear <= 1986 then tr=.46;
if fyear=1987 then tr=.4;
if 1988 <= fyear <= 1992 then tr=.34;
if 1993 <= fyear then tr=.35;
tb_1 = ((txfo + txfed) / tr) / ib;
if missing(txfo) or missing(txfed) then tb_1=((txt-txdi) / tr) / ib; * they rank within industries;
if (txfo + txfed > 0 or txt > txdi) and ib <= 0 then tb_1=1;

txfovars = quarterlycompvars.pivot(columns='LPERMNO', values='txfo')
txfovars = txfovars.fillna(value=None, method='ffill', limit=3)
txfedvars = quarterlycompvars.pivot(columns='LPERMNO', values='txfed')
txfedvars = txfedvars.fillna(value=None, method='ffill', limit=3)
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3 )
txtvars = quarterlycompvars.pivot(columns='LPERMNO', values='txtq')
txtvars = txtvars.fillna(value=None, method='ffill', limit=3)
txdivars = quarterlycompvars.pivot(columns='LPERMNO', values='txdi')
txdivars = txdivars.fillna(value=None, method='ffill', limit=3)
fyearvars = quarterlycompvars.pivot(columns='LPERMNO', values='fyear')
fyearvars = fyearvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
tb_1mask1 = ( ( txfovars.isna() ) | ( txfedvars.isna() ) )
tb_1mask2 = ( ( ( ( txfovars + txfedvars ) > 0 ) | ( txtvars > txdivars ) ) & ( ibvars < 0 ) )

        # Calculating First Variable
tr = pd.DataFrame(columns=txfovars.columns, index=txfovars.index)

trmask1 = ( fyearvars <= 1978)
trmask2 = ( (fyearvars >= 1979) & (fyearvars <= 1986) )
trmask3 = ( fyearvars == 1987 )
trmask4 = ( (fyearvars >= 1988) & (fyearvars <= 1992) )
trmask5 = ( fyearvars >= 1993 )

tr[trmask1] = .48
tr[trmask2] = .46
tr[trmask3] = .44
tr[trmask4] = .34
tr[trmask5] = .35

tb_1 = ( ( txfovars + txfedvars ) / tr.replace(0,np.nan) ) / ibvars.replace(0,np.nan)

         # Applying Mask
tb_1[tb_1mask1] = ( ( txtvars - txdivars ) / tr.replace(0,np.nan) ) / ibvars.replace(0,np.nan)
tb_1[tb_1mask2] = 1

tb_1.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/19.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 20 roa
'''roa = ni / ((at + lag(at)) / 2
For Mohanram score ( 2005 )
'''
nivars = quarterlycompvars.pivot(columns='LPERMNO', values='niq')
nivars = nivars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

roa = nivars / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

roa.astype(float).to_csv(pathtoaddvars_200_400 + '220.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 21 cfroa
'''cfroa = oancf / ((at + lag(at)) / 2);
if missing(oancf) then cfroa=(ib+dp) / ((at+lag(at)) / 2);
For Mohanram score ( 2005 )
'''
oancfvars = quarterlycompvars.pivot(columns='LPERMNO', values='oancfy')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
cfroamask = (oancfvars.isna())

        # Calculating First Variable
cfroa = oancfvars / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

        # Applying Mask
ibvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibvars = ibvars.fillna(value=None, method='ffill', limit=3)
dpvars = quarterlycompvars.pivot(columns='LPERMNO', values='dpq')
dpvars = dpvars.fillna(value=None, method='ffill', limit=3)

cfroa[cfroamask] = ( ibvars + dpvars ) / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

cfroa.astype(float).to_csv(pathtoaddvars_200_400 + '221.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 22 xrdint
'''xrdint = xrd / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )
'''
xrdvars = quarterlycompvars.pivot(columns='LPERMNO', values='xrdq')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

xrdint = xrdvars / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

xrdint.astype(float).to_csv(pathtoaddvars_200_400 + '222.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 23 capxint
'''capxint = capx / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )
'''
capxvars = quarterlycompvars.pivot(columns='LPERMNO', values='capxy')
capxvars = capxvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

capxint = capxvars / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

capxint.astype(float).to_csv(pathtoaddvars_200_400 + '223.csv')

#----------------------------------------------------------------------------------------------------------------------------------
'''XAD IS NOT BEING FOUND IN THE QUARTERLY DATABASE
# Additional Variable 24 xadint
xadint = xad / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )

xadvars = quarterlycompvars.pivot(columns='LPERMNO', values='xad')
xadvars = xadvars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)

xadint = xadvars / (( atvars + atvars.shift(3) ) / 2 ).replace(0,np.nan)

xadint.astype(float).to_csv('/home/emre/Masaüstü/emredenemedosyası/DownGradedAnomalies/24.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 55 orgcap( Organizational Capital )
'''consumer price index to create orgcap measure
from Bureau of Labor Statistics website
retain orgcap_1
avgat=((at+lag(at))/2)
if first.gvkey then orgcap_1=(xsga/cpi)/(.1+.15)
else orgcap_1=orgcap_1*(1-.15)+xsga/cpi
orgcap=orgcap_1/avgat
if count=1 then orgcap=.
'''

cpi = pd.read_csv(pathto_raw_data + 'cpi.csv')
dr = pd.date_range(start= ('1/1/' + str(cpi['Year'][0])), freq='M', periods=len(cpi) )
cpi.index = dr
cpi.index = pd.to_datetime(cpi.index).to_period('M').asfreq('D')
cpi = cpi['Value']

xsga0vars = quarterlycompvars.pivot(columns='LPERMNO', values='xsgaq')
xsga0vars = xsga0vars.fillna(value=None, method='ffill', limit=3)
atvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atvars = atvars.fillna(value=None, method='ffill', limit=3)
countvars = quarterlycompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
orgcapmask1 = countvars==1
orgcapmask2 = countvars>1

        # Calculating First Variable
orgcap_1 = pd.DataFrame(index=countvars.index, columns=countvars.columns)
orgcap_2 = pd.DataFrame(index=countvars.index, columns=countvars.columns)
        # Applying Mask
orgcap_1[orgcapmask1] = ( (xsga0vars.divide(pd.DataFrame(cpi).loc[xsga0vars.index].values, axis=0)) / (.1 + .15) )
orgcap_1 = orgcap_1.fillna(value=None, method='ffill')
orgcap_2[orgcapmask2] =  (xsga0vars.divide(pd.DataFrame(cpi).loc[xsga0vars.index].values, axis=0))
orgcap_3 = orgcap_1 * (1-.15) + orgcap_2
orgcap_3[orgcapmask1] = 0

avgat = ( atvars + atvars.shift(12) ) / 2

orgcap = orgcap_1 / avgat.replace(0,np.nan)

orgcap.to_csv(pathto_200_400 + '255.csv')

#----------------------------------------------------------------------------------------------------------------------------------

#I am Skipping a part which I noted at the Misc File 233 - 241 

#----------------------------------------------------------------------------------------------------------------------------------
    # Description & Function for Industry Adjustments
''' Industry Adjustments for Variables 
						/*other preparation steps for annual variables: industry adjustments*/
						proc sql;
							create table data2
							as select *,chpm-mean(chpm) as chpmia,chato-mean(chato) as chatoia,
							sum(sale) as indsale,hire-mean(hire) as chempia,bm-mean(bm) as bm_ia,
							pchcapx-mean(pchcapx) as pchcapx_ia,tb_1-mean(tb_1) as tb,
							cfp-mean(cfp) as cfp_ia,mve_f-mean(mve_f) as mve_ia
						from data2
						group by sic2,fyear;
						quit;
						proc sql;
						create table data2
						as select *,sum( (sale/indsale)*(sale/indsale) ) as herf
						from data2
						group by sic2,fyear;
						quit;
'''
def indadj(df):

    sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')
    sic = sic.replace(['na', 'Z'], np.nan)
    df_final = pd.DataFrame(index=df.index, columns=df.columns)

    for i in df.index:

        a = df.loc[i].dropna().astype(float)
        b = sic.loc[i].astype(float)
        a.name = 'values'
        b.name = 'sic'
        c = pd.concat([a, b], axis=1)
        d = c.groupby(by='sic').mean()

        c['PERMNO'] = c.index
        c.set_index('sic', inplace=True)
        c = c[c.index.notnull()]
        c = c.sort_index()
        d = d.sort_index()

        c = c[c['values'].notnull()]
        d.dropna(inplace=True)


        c['values2'] = c['values'].astype(float) - d['values'].astype(float)
        c.set_index('PERMNO', inplace=True)

        df_final.loc[i,c.index] = c['values2']

    return df_final

#----------------------------------------------------------------------------------------------------------------------------------
    # 21 chpmia ( Industry-Adjusted Change in Profit Margin )
''' See Description for Industry Adjustments
Use additional Variable 4 ( chpm ) 
'''
chpm = pd.read_csv(pathtoaddvars_200_400 + '204.csv')

chpm.set_index('datadate', inplace=True)
chpm.index = pd.to_datetime(chpm.index).to_period('D')
chpm.columns = chpm.columns.map(int)

chpmia = indadj(chpm)

chpmia.astype(float).to_csv(pathto_200_400 + '221.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 16 chatoia ( Industry-Adjusted Change in Asset Turnover )
''' See Description for Industry Adjustments
Use additional Variable 5 ( chato ) 
'''
chato = pd.read_csv(pathtoaddvars_200_400 + '205.csv')

chato.set_index('datadate', inplace=True)
chato.index = pd.to_datetime(chato.index).to_period('D')
chato.columns = chato.columns.map(int)

chatoia = indadj(chato)

chatoia.astype(float).to_csv(pathto_200_400 + '216.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 25 indsale
''' See Description for Industry Adjustments
'''
salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)

sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')
sic = sic.replace(['na', 'Z'], np.nan)

indsale = pd.DataFrame(index=salevars.index, columns=salevars.columns)

for i in salevars.index:

    a = salevars.loc[i].dropna().astype(float)
    b = sic.loc[i].astype(float)
    a.name = 'values'
    b.name = 'sic'
    c = pd.concat([a, b], axis=1)
    d = c.groupby(by='sic').sum()

    c['PERMNO'] = c.index
    c.set_index('sic', inplace=True)
    c = c[c.index.notnull()]
    c = c.sort_index()
    d = d.sort_index()

    c = c[c['values'].notnull()]
    d.dropna(inplace=True)

    c['values2'] = d['values'].astype(float)
    c.set_index('PERMNO', inplace=True)

    indsale.loc[i, c.index] = c['values2']

indsale.astype(float).to_csv(pathtoaddvars_200_400 + '225.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' # 18 chempia ( Industry-Adjusted Change in Employees )
See Description for Industry Adjustments
Use Variable 38 ( hire )  

hire = pd.read_csv('./data/Anomalies/38.csv')

hire.set_index('datadate', inplace=True)
hire.index = pd.to_datetime(hire.index).to_period('D')
hire.columns = hire.columns.map(int)

chempia = indadj(hire)

chempia.astype(float).to_csv('./data/Anomalies/18.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 10 bm_ia ( Industry-Adjusted Book to Market )
''' See Description for Industry Adjustments
Use Variable 9 ( bm ) 
'''
bm = pd.read_csv(pathto_200_400 + '209.csv')

bm.set_index('datadate', inplace=True)
bm.index = pd.to_datetime(bm.index).to_period('D')
bm.columns = bm.columns.map(int)


bm_ia = indadj(bm)

bm_ia.astype(float).to_csv(pathto_200_400 + '210.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 56 pchcapx_ia ( Industry-Adjusted % chang in Capital Expenditures )
'''See Description for Industry Adjustments
Use Additional Variable 8 ( pchcapx ) 
'''
pchcapx = pd.read_csv(pathtoaddvars_200_400 + '208.csv')

pchcapx.set_index('datadate', inplace=True)
pchcapx.index = pd.to_datetime(pchcapx.index).to_period('D')
pchcapx.columns = pchcapx.columns.map(int)

pchcapx_ia = indadj(pchcapx)

pchcapx_ia.astype(float).to_csv(pathto_200_400 + '256.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' # 92 tb ( Tax Income to Book Income )
See Description for Industry Adjustments
Use Additional Variable 19 ( tb_1 ) 

tb_1 = pd.read_csv('./data/AdditionalVariables/19.csv')

tb_1.set_index('datadate', inplace=True)
tb_1.index = pd.to_datetime(tb_1.index).to_period('D')
tb_1.columns = tb_1.columns.map(int)

tb = indadj(tb_1)

tb.astype(float).to_csv('./data/Anomalies/92.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    # 15 cfp_ia ( Industry_Adjusted Cash Flow to Price Ratio )
'''See Description for Industry Adjustments
Use Variable 14 ( cfp ) 
'''
cfp = pd.read_csv(pathto_200_400 + '214.csv')

cfp.set_index('datadate', inplace=True)
cfp.index = pd.to_datetime(cfp.index).to_period('D')
cfp.columns = cfp.columns.map(int)

cfp_ia = indadj(cfp)

cfp_ia.astype(float).to_csv(pathto_200_400 + '215.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' SAME VARIABLE 
 # 52 mve_ia ( Industry_Adjusted size )
See Description for Industry Adjustments

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

mve_ia = indadj(mve_f)

mve_ia.astype(float).to_csv('./data/Anomalies/52.csv')
'''
#----------------------------------------------------------------------------------------------------------------------------------
    #  37 herf ( Industry Sale Concentration  ) 
'''See Description for Industry Adjustments
sum( (sale/indsale)*(sale/indsale) ) as herf
'''
indsale = pd.read_csv(pathtoaddvars_200_400 + '225.csv')

indsale.set_index('datadate', inplace=True)
indsale.index = pd.to_datetime(indsale.index).to_period('D')
indsale.columns = indsale.columns.map(int)

salevars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
salevars = salevars.fillna(value=None, method='ffill', limit=3)

sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')
sic = sic.replace(['na', 'Z'], np.nan)

herf = ( salevars / indsale.replace(0,np.nan) ) ** 2

for i in herf.index:

    a = herf.loc[i].dropna().astype(float)
    b = sic.loc[i].astype(float)
    a.name = 'values'
    b.name = 'sic'
    c = pd.concat([a, b], axis=1)
    d = c.groupby(by='sic').sum()

    c['PERMNO'] = c.index
    c.set_index('sic', inplace=True)
    c = c[c.index.notnull()]
    c = c.sort_index()
    d = d.sort_index()

    c = c[c['values'].notnull()]
    d.dropna(inplace=True)

    c['values2'] = d['values'].astype(float)
    c.set_index('PERMNO', inplace=True)

    herf.loc[i, c.index] = c['values2']

herf.astype(float).to_csv(pathto_200_400 + '237.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 50 Financial Statement Score
'''Mohanram Score
ms = m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8
if roa > md_roa then m1=1; else m1=0
if cfroa > md_cfroa then m2=1; else m2=0
if oancf > ni then m3=1; else m3=0
if xrdint > md_xrdint then m4=1; else m4=0
if capxint > md_capxint then m5=1; else m5=0
if xadint > xadint then m6=1; else m6=0
'''
def ia_md(df):

    sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')
    sic = sic.replace(['na', 'Z'], np.nan)
    
    final_df = pd.DataFrame(index=df.index, columns=df.columns)
    
    for i in df.index:
    
        a = df.loc[i].dropna().astype(float)
        b = sic.loc[i].astype(float)
        a.name = 'values'
        b.name = 'sic'
        c = pd.concat([a, b], axis=1)
        d = c.groupby(by='sic').median()
    
        c['PERMNO'] = c.index
        c.set_index('sic', inplace=True)
        c = c[c.index.notnull()]
        c = c.sort_index()
        d = d.sort_index()
    
        c = c[c['values'].notnull()]
        d.dropna(inplace=True)
    
        c['values2'] = d['values'].astype(float)
        c.set_index('PERMNO', inplace=True)
    
        final_df.loc[i, c.index] = c['values2']

    return final_df

    # M1
roa = pd.read_csv(pathtoaddvars_1_200 + '20.csv')

roa.set_index('datadate', inplace=True)
roa.index = pd.to_datetime(roa.index).to_period('D')
roa.columns = roa.columns.map(int)

#md_roa = ia_md(roa)
#md_roa.to_csv(pathtoaddvars_1_200 + 'm1.csv')
md_roa = pd.read_csv(pathtoaddvars_1_200 + 'm1.csv')
md_roa.set_index('datadate', inplace=True)
md_roa.index = pd.to_datetime(md_roa.index).to_period('D')
md_roa.columns = md_roa.columns.map(int)

m1mask = ( roa > md_roa )
m1mask2 = ( roa < md_roa )

m1 = pd.DataFrame(index=roa.index, columns=roa.columns )

m1[m1mask] = 1 
m1[m1mask2] = 0

    # M2
cfroa = pd.read_csv(pathtoaddvars_1_200 + '21.csv')

cfroa.set_index('datadate', inplace=True)
cfroa.index = pd.to_datetime(cfroa.index).to_period('D')
cfroa.columns = cfroa.columns.map(int)

#md_cfroa = ia_md(cfroa)    
#md_cfroa.to_csv(pathtoaddvars_1_200 + 'm2.csv')
md_cfroa = pd.read_csv(pathtoaddvars_1_200 + 'm2.csv')
md_cfroa.set_index('datadate', inplace=True)
md_cfroa.index = pd.to_datetime(md_cfroa.index).to_period('D')
md_cfroa.columns = md_cfroa.columns.map(int)    

m2mask = ( cfroa > md_cfroa )
m2mask2 = ( cfroa < md_cfroa )

m2 = pd.DataFrame(index=cfroa.index, columns=cfroa.columns )

m2[m2mask] = 1   
m2[m2mask2] = 0

    # M3
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
nivars = annualcompvars.pivot(columns='LPERMNO', values='ni')
nivars = nivars.fillna(value=None, method='ffill', limit=12)

m3mask = ( oancfvars > nivars )
m3mask2 = ( oancfvars < nivars )

m3 = pd.DataFrame(index=oancfvars.index, columns=oancfvars.columns )

m3[m3mask] = 1 
m3[m3mask2] = 0

    # M4
xrdint = pd.read_csv(pathtoaddvars_1_200 + '22.csv')

xrdint.set_index('datadate', inplace=True)
xrdint.index = pd.to_datetime(xrdint.index).to_period('D')
xrdint.columns = xrdint.columns.map(int)

#md_xrdint = ia_md(xrdint)    
#md_xrdint.to_csv(pathtoaddvars_1_200 + 'm4.csv')
md_xrdint = pd.read_csv(pathtoaddvars_1_200 + 'm4.csv')
md_xrdint.set_index('datadate', inplace=True)
md_xrdint.index = pd.to_datetime(md_xrdint.index).to_period('D')
md_xrdint.columns = md_xrdint.columns.map(int)    

m4mask = ( xrdint > md_xrdint )
m4mask2 = ( xrdint < md_xrdint )

m4 = pd.DataFrame(index=xrdint.index, columns=xrdint.columns )

m4[m4mask] = 1 
m4[m4mask2] = 0

    # M5
capxint = pd.read_csv(pathtoaddvars_1_200 + '23.csv')

capxint.set_index('datadate', inplace=True)
capxint.index = pd.to_datetime(capxint.index).to_period('D')
capxint.columns = capxint.columns.map(int)

#md_capxint = ia_md(capxint)    
#md_capxint.to_csv(pathtoaddvars_1_200 + 'm5.csv')
md_capxint = pd.read_csv(pathtoaddvars_1_200 + 'm5.csv')
md_capxint.set_index('datadate', inplace=True)
md_capxint.index = pd.to_datetime(md_capxint.index).to_period('D')
md_capxint.columns = md_capxint.columns.map(int)  

m5mask = ( capxint > md_capxint )
m5mask2 = ( capxint < md_capxint )

m5 = pd.DataFrame(index=capxint.index, columns=capxint.columns )

m5[m5mask] = 1 
m5[m5mask2] = 0

    # M6
xadint = pd.read_csv(pathtoaddvars_1_200 + '24.csv')

xadint.set_index('datadate', inplace=True)
xadint.index = pd.to_datetime(xadint.index).to_period('D')
xadint.columns = xadint.columns.map(int)

#md_xadint = ia_md(xadint)    
#md_xadint.to_csv(pathtoaddvars_1_200 + 'm6.csv')
md_xadint = pd.read_csv(pathtoaddvars_1_200 + 'm6.csv')
md_xadint.set_index('datadate', inplace=True)
md_xadint.index = pd.to_datetime(md_xadint.index).to_period('D')
md_xadint.columns = md_xadint.columns.map(int) 

m6mask = ( xadint > md_xadint )
m6mask2 = ( xadint < md_xadint )

m6 = pd.DataFrame(index=xadint.index, columns=xadint.columns )

m6[m6mask] = 1 
m6[m6mask2] = 0

    # M7 
roavol = pd.read_csv(pathto_1_200 + '75.csv')

roavol.set_index('datadate', inplace=True)
roavol.index = pd.to_datetime(roavol.index).to_period('D')
roavol.columns = roavol.columns.map(int)

#md_roavol = ia_md(roavol) 
#md_roavol.to_csv(pathtoaddvars_1_200 + 'm7.csv')

md_roavol = pd.read_csv(pathtoaddvars_1_200 + 'm7.csv')
md_roavol.set_index('datadate', inplace=True)
md_roavol.index = pd.to_datetime(md_roavol.index).to_period('D')
md_roavol.columns = md_roavol.columns.map(int) 

m7mask = ( roavol > md_roavol )
m7mask2 = ( roavol < md_roavol )

m7 = pd.DataFrame(index=roavol.index, columns=roavol.columns )

m7[m7mask] = 1 
m7[m7mask2] = 0 

    # M8 
sgrvol = pd.read_csv(pathtoaddvars_1_200 + '29.csv')

sgrvol.set_index('datadate', inplace=True)
sgrvol.index = pd.to_datetime(sgrvol.index).to_period('D')
sgrvol.columns = sgrvol.columns.map(int)

#md_sgrvol = ia_md(sgrvol) 
#md_sgrvol.to_csv(pathtoaddvars_1_200 + 'm8.csv')
md_sgrvol = pd.read_csv(pathtoaddvars_1_200 + 'm8.csv')
md_sgrvol.set_index('datadate', inplace=True)
md_sgrvol.index = pd.to_datetime(md_sgrvol.index).to_period('D')
md_sgrvol.columns = md_sgrvol.columns.map(int) 


m8mask = ( sgrvol > md_sgrvol )
m8mask2 = ( sgrvol < md_sgrvol )

m8 = pd.DataFrame(index=sgrvol.index, columns=sgrvol.columns )

m8[m8mask] = 1 
m8[m8mask2] = 0

ms = m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8

countvars = monthlycrspvars.pivot(columns='PERMNO', values='count')
countvars[countvars.notnull()] = 1

ms = ms*countvars
ms.astype(float).to_csv(pathto_200_400 + '250.csv')
'''

'''*this is for abnormal trading volume and returns around earings announcements;
proc sql;
	create table data5 
	as select a.*,b.vol
	from data5 a left join crsp.dsf b
	on a.permno=b.permno and
     intnx('WEEKDAY',rdq,-30)<=b.date<=intnx('WEEKDAY',rdq,-10);
	quit; 	
						proc sql;
							create table data5
							as select *,mean(vol) as avgvol
							from data5
						group by permno,datadate,rdq;
						quit;
						proc sort data=data5(drop=vol) nodupkey;
						where not missing(rdq);
						by permno datadate rdq;
						run;    									
						proc sql;
						create table data6 
						as select a.*,b.vol,b.ret
						from data5 a left join crsp.dsf b
						on a.permno=b.permno and
   						  intnx('WEEKDAY',rdq,-1)<=b.date<=intnx('WEEKDAY',rdq,1);
						quit;
						proc sql;
						create table data6
						as select *,(mean(vol)-avgvol)/avgvol as aeavol,sum(ret) as ear
						from data6
						group by permno,datadate,rdq;
						quit;
						proc sort data=data6(drop=vol avgvol ret) nodupkey;
						by permno datadate rdq;
						run;
'''

#----------------------------------------------------------------------------------------------------------------------------------
    # 3 aeavol ( Abnormal Earnings Announcement Volume )
rdqvars = quarterlycompvars[['LPERMNO', 'rdq']]
rdqvars.reset_index(inplace=True)
rdqvars = rdqvars.loc[rdqvars['rdq'].dropna().index]
rdqvars['dum'] = 1
rdq = rdqvars.drop_duplicates(subset=['rdq','LPERMNO']).pivot(index='rdq', columns='LPERMNO', values='dum')
rdq.index = pd.to_datetime(rdq.index, format='%d/%m/%Y')
rdq = rdq.sort_index()
rdq = rdq.loc[rdq.index.dropna()]

volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

volmean20 = volvars.rolling(20).mean().shift(10)
volmean3 = volvars.rolling(3).mean().shift(-1)

avgvolvars20 = volmean20.loc[volmean20.index[(rdq.index[:-109]).map(lambda x: volmean20.index.get_loc(x, method='bfill'))]]
avgvolvars3 = volmean3.loc[volmean3.index[(rdq.index[:-109]).map(lambda x: volmean3.index.get_loc(x, method='bfill'))]]

rdq2_1 = rdq.iloc[:-109].copy() 
rdq2_1.index = volmean20.index[(rdq.index[:-109]).map(lambda x: volmean20.index.get_loc(x, method='bfill'))]

rdq2_2 = rdq.iloc[:-109].copy() 
rdq2_2.index = volmean3.index[(rdq.index[:-109]).map(lambda x: volmean3.index.get_loc(x, method='bfill'))]

avgvol20 = rdq2_1 * avgvolvars20.loc[:,rdq2_1.columns]
avgvol3 = rdq2_2 * avgvolvars3.loc[:,rdq2_2.columns]

aeavol = (avgvol3 - avgvol20 ) /avgvol20.replace(0,np.nan)
aeavol.index = pd.to_datetime(aeavol.index).to_period('M').asfreq('D')
aeavol = aeavol.groupby(by=aeavol.index).sum()
aeavol = aeavol.replace(0,np.nan)

aeavol.index = (aeavol.index).to_timestamp().to_period('M')
aeavol = aeavol.groupby(by=aeavol.index).last()
aeavol.index = (aeavol.index).asfreq('D')

aeavol.to_csv(pathto_200_400 + '203.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 31 ear ( Earnings Announcement Return )
rdqvars = quarterlycompvars[['LPERMNO', 'rdq']]
rdqvars.reset_index(inplace=True)
rdqvars = rdqvars.loc[rdqvars['rdq'].dropna().index]
rdqvars['dum'] = 1
rdq = rdqvars.drop_duplicates(subset=['rdq','LPERMNO']).pivot(index='rdq', columns='LPERMNO', values='dum')
rdq.index = pd.to_datetime(rdq.index, format='%d/%m/%Y')
rdq = rdq.sort_index()
rdq = rdq.loc[rdq.index.dropna()]

retvars = pd.read_csv('/run/media/research/Kaan/memoryerrortransferringfile/WRDS_Data_Preprocessing/data/WRDSdata/092019/13092019/dailyret.csv')
retvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
retvars = retvars.pivot(index='date', columns='PERMNO', values='RET')
retvars.index = pd.to_datetime(retvars.index, format='%d/%m/%Y')
retvars = retvars.sort_index()

retvars = retvars.replace('C', np.nan)
retvars = retvars.replace('B', np.nan)
retvars = retvars.astype('float')

cumret = (retvars+1).rolling(3).apply(np.prod)-1

earvars = cumret.loc[cumret.index[(rdq.index[:-109]).map(lambda x: cumret.index.get_loc(x, method='bfill'))]]

rdq2 = rdq.iloc[:-109].copy() 
rdq2.index = cumret.index[(rdq.index[:-109]).map(lambda x: cumret.index.get_loc(x, method='bfill'))]

ear = rdq2 * earvars

ear.index = pd.to_datetime(ear.index).to_period('M').asfreq('D')
ear = ear.replace(np.nan, 0)
ear = ear.groupby(by=ear.index).sum()

ear.index = (ear.index).to_timestamp().to_period('M')
ear = ear.groupby(by=ear.index).last()
ear.index = (ear.index).asfreq('D')

ear = ear.replace(0, np.nan)

ear.to_csv(pathto_200_400 + '231.csv')
