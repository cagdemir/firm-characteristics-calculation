import pandas as pd
import numpy as np
import multiprocessing as mp

pathtorawdata = '...'

#----------------------------------------------------------------------------------------------------------------------------------

pathto_1_200 = '...'
pathto_200_400 = '...'
pathto_400_500 = '...'
pathto_500_600 = '...'

#----------------------------------------------------------------------------------------------------------------------------------

pathtoaddvars_1_200 = '...'
pathtoaddvars_200_400 = '...'

#----------------------------------------------------------------------------------------------------------------------------------
''' COMPUSTAT ANNUAL INFORMATION '''

### Firm Variables

    # Header Information
'''substr(compress(cusip),1,6) as cnum,c.gvkey,datadate,fyear,c.cik,substr(sic,1,2) as sic2,sic,naics'''
'''Cusip: A cusip number is a 9-character alphanumeric code that identifies a security for the purposes of facilitating clearing and settlement of trades.   
       Gvkey: Standart and Poor's Identifier.
       Datadate:
       Fyear: Fiscal Year ( year -1 if the date falls between 01/01 - 01/06, year if the date falls between 01/06 - 01/01. )
       cik: Central Index key ( it is number given to an individual, company, of foreign government by the US Securities and Exchange Commission. )
       sic: Standard Industry Classification Code.
       naics: North American Industry Classification Code.       
       '''

    # Income Statement
'''sale,revt,cogs,xsga,dp,xrd,xad,ib,ebitda,ebit,nopi,spi,pi,txp,ni,txfed,txfo,txt,xint
       sale: The amount of actual billings to customers for regular sales completed during the period.
       revt: Total Revenue
       cogs: Cost of Goods Sold
       xsga: Selling, General and Administrative Expense
       dp: Depreciation and Amortization Expense
       xrd: Research and Development Expense
       xad: Advertising Expense
       ib: Income Before Extraordinary Items
       ebitda:  Earnings Before Interest
       ebit: Earnings Before Interest and Taxes
       nopi: Nonoperating Income (Expense)
       spi: Special Items
       pi: Pretax Income
       txp: Income Taxes Payable
       ni: Net Income (Loss)
       txfed: Income Taxes Federal
       txfo: Income Taxes - Foreign
       txt: Income Taxes - Total
       xint: Interest and Related Expense - Total
    '''

    # CF statement and others
'''capx,oancf,dvt,ob,gdwlia,gdwlip,gwo'''

    # Assets
'''rect,act,che,ppegt,invt,at,aco,intan,ao,ppent,gdwl,fatb,fatl'''

    # Liabilities
'''lct,dlc,dltt,lt,dm,dcvt,cshrc,dcpstk,pstk,ap,lco,lo,drc,drlt,txdi'''

    # Equity and other
'''ceq,scstkc,emp,csho
       ceq: Total common/ordinary equity.
    '''

    # Market
'''abs(prcc_f) as prcc_f,csho*calculated prcc_f as mve_f'''


#----------------------------------------------------------------------------------------------------------------------------------

annualcompvars = pd.read_csv((pathtorawdata + 'annualcomp.csv'), usecols=['LPERMNO','cusip','GVKEY','datadate','fyear','cik','sic','naics','sale','revt','cogs','xsga','dp','xrd','xad','ib','ebitda','ebit','nopi','spi','pi','txp','ni','txfed','txfo','txt','xint','capx','oancf','dvt','ob','gdwlia','gdwlip','gwo','rect','act','che','ppegt','invt','at','aco','intan','ao','ppent','gdwl','fatb','fatl','lct','dlc','dltt','lt','dm','dcvt','cshrc','dcpstk','pstk','ap','lco','lo','drc','drlt','txdi','ceq','scstkc','emp','csho','prcc_f'])

annualcompvars.drop_duplicates(subset=['datadate', 'LPERMNO'], inplace=True)
annualcompvars.set_index('datadate', inplace=True)
annualcompvars.index = pd.to_datetime(annualcompvars.index, format='%d/%m/%Y').to_period('M').asfreq('D')

#----------------------------------------------------------------------------------------------------------------------------------

''' COMPUSTAT QUARTERLY INFORMATION'''
    ### Firm Variables

    # Header Information
'''select substr(compress(cusip),1,6) as cnum,c.gvkey,fyearq,fqtr,datadate,rdq,substr(sic,1,2) as sic2'''

    # Income Statement
'''ibq,saleq,txtq,revtq,cogsq,xsgaq'''

    # Balance Sheet Items
'''atq,actq,cheq,lctq,dlcq,ppentq'''

    # Others
'''abs(prccq) as prccq,abs(prccq)*cshoq as mveq,ceqq, seqq,pstkq,atq,ltq,pstkrq'''

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

quarterlycompvars = pd.read_csv((pathtorawdata + 'quartcomp.csv'), usecols=['LPERMNO','cusip','GVKEY','datadate','cik','sic','naics','ibq', 'rdq','saleq','txtq','revtq','cogsq','prccq','cshoq','xsgaq','atq','actq','cheq','lctq','dlcq','ppentq','ceqq','seqq','pstkq','atq','ltq','pstkrq'])

quarterlycompvars.drop_duplicates(subset=['datadate', 'LPERMNO'], inplace=True)
quarterlycompvars.set_index('datadate', inplace=True)
quarterlycompvars.index = pd.to_datetime(quarterlycompvars.index, format='%d/%m/%Y').to_period('M').asfreq('D')

#----------------------------------------------------------------------------------------------------------------------------------

'''CRSP MONTHLY INFORMATION'''

monthlycrspvars = pd.read_csv((pathtorawdata + 'monthlycrsp.csv'), usecols=['PERMNO','date', 'SHROUT', 'ALTPRC', 'CUSIP', 'EXCHCD', 'SHRCD', 'SICCD', 'SPREAD','VOL'])

monthlycrspvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
monthlycrspvars.set_index('date', inplace=True)
monthlycrspvars.index = pd.to_datetime(monthlycrspvars.index, format='%d/%m/%Y').to_period('M').asfreq('D')

#----------------------------------------------------------------------------------------------------------------------------------

    # Market Cap

        # Calculating Annual Market Cap

annualcompvars['mve_f'] = annualcompvars['prcc_f'].abs() * annualcompvars['csho']
annualcompvars['mve_f'].replace(0, np.nan, inplace=True)

        # Calculating Quarterly Market Cap

quarterlycompvars['mve_f'] = quarterlycompvars['prccq'].abs() * quarterlycompvars['cshoq']
quarterlycompvars['mve_f'].replace(0, np.nan, inplace=True)


        # Calculating Monthly Market Cap

monthlycrspvars['mve_f'] = monthlycrspvars['ALTPRC'].abs() * monthlycrspvars['SHROUT']
monthlycrspvars['mve_f'].replace(0, np.nan, inplace=True)

#----------------------------------------------------------------------------------------------------------------------------------

    # Compress cusip 1:6
annualcompvars['cnum'] = annualcompvars['cusip'].apply(lambda x: str(x)[:8])
quarterlycompvars['cnum'] = quarterlycompvars['cusip'].apply(lambda x: str(x)[:8])

    # Substr sic 1:2
annualcompvars['sic2'] = annualcompvars['sic'].apply(lambda x: str(x)[:2])
quarterlycompvars['sic2'] = quarterlycompvars['sic'].apply(lambda x: str(x)[:2])
monthlycrspvars['SICCD2'] = monthlycrspvars['SICCD'].apply(lambda x: str(x)[:2])

#----------------------------------------------------------------------------------------------------------------------------------
    # Count
annualcompvars['count'] = 1
quarterlycompvars['count'] = 1
monthlycrspvars['count'] = 1

annualcompvars['count'] = pd.DataFrame(annualcompvars.groupby(by='GVKEY').cumsum()['count'])
quarterlycompvars['count'] = pd.DataFrame(quarterlycompvars.groupby(by='GVKEY').cumsum()['count'])
monthlycrspvars['count'] = pd.DataFrame(monthlycrspvars.groupby(by='PERMNO').cumsum()['count'])



'''						*do some clean up, several of these variables have lots of missing values;
						if not missing(drc) and not missing(drlt) then dr=drc+drlt;
						if not missing(drc) and missing(drlt) then dr=drc;
						if not missing(drlt) and missing(drc) then dr=drlt;
						
						if missing(dcvt) and not missing(dcpstk) and not missing(pstk) and dcpstk>pstk then dc=dcpstk-pstk;
						if missing(dcvt) and not missing(dcpstk) and missing(pstk) then dc=dcpstk;
						if missing(dc) then dc=dcvt;

						if missing(xint) then xint0=0;
							else xint0=xint;
						if missing(xsga) then xsga0=0;
							else xsga0=0;
						
							'''

    # dr calculation
annualcompvars['dr'] = annualcompvars['drc'].replace(np.nan, 0) + annualcompvars['drlt'].replace(np.nan, 0)

    # dc calculation
annualcompvars['dc'] = np.nan

dc_mask_1 = (annualcompvars['dcvt'].isna() & annualcompvars['dcpstk'].notnull() & annualcompvars['pstk'].notnull() & annualcompvars['dcpstk']>annualcompvars['pstk'])
annualcompvars['dc'][dc_mask_1] = annualcompvars['dcpstk'] - annualcompvars['pstk']

dc_mask_2 = (annualcompvars['dcvt'].isna() & annualcompvars['dcpstk'].notnull() & annualcompvars['pstk'].isna())
annualcompvars['dc'][dc_mask_2] = annualcompvars['dcpstk']

dc_mask_3 = (annualcompvars['dc'].isna())
annualcompvars['dc'][dc_mask_3] = annualcompvars['dcvt']

    # xint0 calculation
annualcompvars['xint0'] = annualcompvars['xint'].replace(np.nan,0)

    # xsga0 calculation
annualcompvars['xsga0'] = annualcompvars['xsga'].replace(np.nan,0)



#----------------------------------------------------------------------------------------------------------------------------------

### create simple-just annual Compustat variables

#----------------------------------------------------------------------------------------------------------------------------------
    # 9 bm ( Book to Market )
'''book-to-Market = (Common Equity / Market Cap)
'''
bmvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
bmvars = bmvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

bm = bmvars / mve_f

bm.astype(float).to_csv(pathto_1_200 + '9.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 33 ep ( Earnings to Price )
'''Earning-to-Price = ( Income Before Extra Ordinary Items / Market Cap )
'''
epvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
epvars = epvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

ep = epvars / mve_f

ep.astype(float).to_csv(pathto_1_200 + '33.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 13 cashpr ( Cash Productivity )
'''Cash Productivity = ( ( Market Cap + Long Term Debt - Total Asset) / Cash and Short Term Investments)
'''
ltdebtvars = annualcompvars.pivot(columns='LPERMNO', values='dltt')
ltdebtvars = ltdebtvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
chevars = chevars.replace(0,np.nan)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

cashpr = ((mve_f + ltdebtvars - atvars) / chevars)

cashpr.astype(float).to_csv(pathto_1_200 + '13.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 30 dy (Dividend to Price)
'''dy = dvt / mve_f
'''
dvtvars = annualcompvars.pivot(columns='LPERMNO', values='dvt')
dvtvars = dvtvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

dy = dvtvars / mve_f

dy.astype(float).to_csv(pathto_1_200 + '30.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 43 lev ( Leverage )
'''lev = lt / mve_f
'''
levvars = annualcompvars.pivot(columns='LPERMNO', values='lt')
levvars = levvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

lev = levvars / mve_f

lev.astype(float).to_csv(pathto_1_200 + '43.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 86 sp ( Sales to Price )
'''sp = sale / mve_f
'''
spvars = annualcompvars.pivot(columns='LPERMNO', values='sale')
spvars = spvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

sp = spvars / mve_f

sp.astype(float).to_csv(pathto_1_200 + '86.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 77 roic ( Return on Ivested Capital )
'''roic = (ebit - nopi) / (ceq + lt - che)
'''
ebitvars = annualcompvars.pivot(columns='LPERMNO', values='ebit')
ebitvars = ebitvars.fillna(value=None, method='ffill', limit=12)
nopivars = annualcompvars.pivot(columns='LPERMNO', values='nopi')
nopivars = nopivars.fillna(value=None, method='ffill', limit=12)
ceqvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=12)
ltvars = annualcompvars.pivot(columns='LPERMNO', values='lt')
ltvars = ltvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)

denominator = ceqvars + ltvars - chevars
denominator = denominator.replace(0,np.nan)

roic = (ebitvars - nopivars) / denominator

roic.astype(float).to_csv(pathto_1_200 + '77.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 71 rd_sale ( R&D to sale )
'''rd_sale = xrd / sale
'''
xrdvars = annualcompvars.pivot(columns='LPERMNO', values='xrd')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=12)
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)

rd_sale = xrdvars / salevars.replace(0,np.nan)

rd_sale.astype(float).to_csv(pathto_1_200 + '71.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 70 rd_mve ( R&D to Market Capitalization )
'''rd_mve = xrd / mve_f
'''
xrdvars = annualcompvars.pivot(columns='LPERMNO', values='xrd')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

rd_mve = xrdvars / mve_f

rd_mve.astype(float).to_csv(pathto_1_200 + '70.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 5 agr ( Asset Growth )
'''agr = (at / lag(at)) - 1
'''
agrvars = annualcompvars.pivot(columns='LPERMNO', values='at')
agrvars = agrvars.fillna(value=None, method='ffill', limit=12)

agr = (agrvars / agrvars.replace(0,np.nan).shift(12)) - 1

agr.astype(float).to_csv(pathto_1_200 + '5.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 34 gma ( Gross Profitability )
'''gma = (revt - cogs) / lag(at)
'''
revtvars = annualcompvars.pivot(columns='LPERMNO', values='revt')
revtvars = revtvars.fillna(value=None, method='ffill', limit=12)
cogsvars = annualcompvars.pivot(columns='LPERMNO', values='cogs')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

gma = (revtvars - cogsvars) / atvars.replace(0,np.nan).shift(12)

gma.astype(float).to_csv(pathto_1_200 + '34.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 17 chcsho ( Change in Shares Outstanding )
'''chcsho = (csho / lag(csho)) - 1
'''
csho = monthlycrspvars.pivot(columns='PERMNO', values='SHROUT')

chcsho = (csho / csho.replace(0,np.nan).shift()) - 1

chcsho.astype(float).to_csv(pathto_1_200 + '17.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 44 lgr ( Growth in Long Term Debt )
'''lgr = (lt / lag(lt)) - 1'''
lgrvars = annualcompvars.pivot(columns='LPERMNO', values='lt')
lgrvars = lgrvars.fillna(value=None, method='ffill', limit=12)

lgr = (lgrvars / lgrvars.replace(0,np.nan).shift(12)) - 1

lgr.astype(float).to_csv(pathto_1_200 + '44.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 2 acc ( Accruals )
'''	acc=(ib-oancf) /  ((at+lag(at))/2);
if missing(oancf) then acc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/  ((at+lag(at))/2)
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

        # Creating mask
accmask = (oancfvars.isna())

        # Calculating First Variable
acc = (ibvars - oancfvars) / ((atvars + atvars.shift(12)).replace(0,np.nan) / 2)

        # Applying mask
actvars = annualcompvars.pivot(columns='LPERMNO', values='act')
actvars = actvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
lctvars = annualcompvars.pivot(columns='LPERMNO', values='lct')
lctvars = lctvars.fillna(value=None, method='ffill', limit=12)
dlcvars = annualcompvars.pivot(columns='LPERMNO', values='dlc')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=12)
txpvars = annualcompvars.pivot(columns='LPERMNO', values='txp')
txpvars = txpvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

acc[accmask] =( ( ( (actvars - actvars.shift(12)) -
                             (chevars - chevars.shift(12)) ) -
                           ( (lctvars - lctvars.shift(12)) -
                             (dlcvars - dlcvars.shift(12)) -
                             (txpvars - txpvars.shift(12)) -
                             (dpvars) ) ) /
                           ( (atvars + atvars.shift(12)).replace(0,np.nan) / 2) )

acc.astype(float).to_csv(pathto_1_200 + '2.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 65 pctacc ( Percent Accruals )
'''pctacc=(ib-oancf)/abs(ib);
if ib=0 then pctacc=(ib-oancf)/.01;
if missing(oancf) then pctacc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/abs(ib);
if missing(oancf) and ib=0 then pctacc=(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) )/.01;
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)

        # Creating Masks
pctaccmask1 = (ibvars==0)
pctaccmask2 = (oancfvars.isna())
pctaccmask3 = ((oancfvars.isna()) & (ibvars==0))

        # Calculating first variable
pctacc = (ibvars - oancfvars) / ibvars.abs()

        # Applying First Mask
pctacc[pctaccmask1] = (ibvars - oancfvars) / .01

        # Applying Second Mask
actvars = annualcompvars.pivot(columns='LPERMNO', values='act')
actvars = actvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
lctvars = annualcompvars.pivot(columns='LPERMNO', values='lct')
lctvars = lctvars.fillna(value=None, method='ffill', limit=12)
dlcvars = annualcompvars.pivot(columns='LPERMNO', values='dlc')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=12)
txpvars = annualcompvars.pivot(columns='LPERMNO', values='txp')
txpvars = txpvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)

pctacc[pctaccmask2] = ( ( ( (actvars - actvars.shift(12)) -
                            (chevars - chevars.shift(12)) ) -
                          ( (lctvars - lctvars.shift(12)) -
                            (dlcvars - dlcvars.shift(12)) -
                            (txpvars - txpvars.shift(12)) -
                            (dpvars) ) ) /
                          ( (ibvars.abs()) ) )

        # Applying Third Mask
pctacc[pctaccmask3] = ( ( ( (actvars - actvars.shift(12)) -
                            (chevars - chevars.shift(12)) ) -
                          ( (lctvars - lctvars.shift(12)) -
                            (dlcvars - dlcvars.shift(12)) -
                            (txpvars - txpvars.shift(12)) -
                            (dpvars) ) ) /
                          ( (.01) ) )

pctacc.astype(float).to_csv(pathto_1_200 + '65.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 14 cfp ( Cash Flow to Price Ratio )
'''cfp=oancf/mve_f
if missing(oancf) then cfp=(ib-(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) ))/mve_f
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
actvars = annualcompvars.pivot(columns='LPERMNO', values='act')
actvars = actvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
lctvars = annualcompvars.pivot(columns='LPERMNO', values='lct')
lctvars = lctvars.fillna(value=None, method='ffill', limit=12)
dlcvars = annualcompvars.pivot(columns='LPERMNO', values='dlc')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=12)
txpvars = annualcompvars.pivot(columns='LPERMNO', values='txp')
txpvars = txpvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

        # Creating mask
cfpmask = (oancfvars.isna())

        # Calculating First Variable
cfp = oancfvars / mve_f

        # Applying mask
cfp[cfpmask] = ( (ibvars -
    ( ( (actvars - actvars.shift(12)) -
        (chevars - chevars.shift(12)) ) -
      ( (lctvars - lctvars.shift(12)) -
        (dlcvars - dlcvars.shift(12)) -
        (txpvars - txpvars.shift(12)) -
        (dpvars) ) ) ) /
         mve_f )


cfp.astype(float).to_csv(pathto_1_200 + '14.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 1 absacc ( Absolute Accruals )
'''Absolute of Accruals
Use Variable #2 acc ( Accruals)
'''
acc = pd.read_csv('./data/Anomalies/2.csv')
acc.set_index('datadate', inplace=True)
acc.index = pd.to_datetime(acc.index).to_period('M').asfreq('D')

absacc = acc.abs()

absacc.astype(float).to_csv(pathto_1_200 + '1.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 4 age ( # of years Since First Compustat Coverage )
'''	if first.gvkey then count=1;
	else count+1;
Original Paper's Description: Firm Age is defined as the number of months between event month t and the first month that a stock appears in CRSP'''
agevars = annualcompvars.pivot(columns='LPERMNO', values='count')
agevars = agevars.fillna(value=None, method='ffill', limit=12)

age = agevars.copy()

age.astype(float).to_csv(pathto_1_200 + '4.csv')
#----------------------------------------------------------------------------------------------------------------------------------
	# 19 chinv ( Change in Inventory )
'''chinv=(invt-lag(invt))/((at+lag(at))/2)
'''
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

chinv = ( invtvars - invtvars.shift(12) ) / ( ( atvars + atvars.shift(12)).replace(0,np.nan) / 2 )

chinv.astype(float).to_csv(pathto_1_200 + '19.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 1 spii
'''if spi ne 0 and not missing(spi) then spii=1; else spii=0
'''
annualcompvars['spii'] = 0
annualcompvars['spii'][(annualcompvars['spi']!=0) & (annualcompvars['spi'].notnull())] = 1

spii = annualcompvars.pivot(columns='LPERMNO', values='spii')
spii = spii.fillna(value=None, method='ffill', limit=12)

spii.astype(float).to_csv(pathtoaddvars_1_200 + '1.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 2 spi ( Altered spi )
'''spi=spi/ ((at+lag(at))/2)
'''
spivars = annualcompvars.pivot(columns='LPERMNO', values='spi')
spivars = spivars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

spi = spivars / ( (atvars + atvars.shift(12)) / 2 ).replace(0,np.nan)

spi.astype(float).to_csv(pathtoaddvars_1_200 + '2.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 3 cf
'''cf=oancf/((at+lag(at))/2);
if missing(oancf) then cf=(ib-(	(act-lag(act) - (che-lag(che))) - (  (lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp ) ))/((at+lag(at))/2)
'''
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

        # Creating mask
cfmask = (oancfvars.isna())

        # Calculating First Variable
cf = oancfvars / ( (atvars + atvars.shift(12)) / 2 ).replace(0,np.nan)

        # Applying Mask
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
actvars = annualcompvars.pivot(columns='LPERMNO', values='act')
actvars = actvars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
lctvars = annualcompvars.pivot(columns='LPERMNO', values='lct')
lctvars = lctvars.fillna(value=None, method='ffill', limit=12)
dlcvars = annualcompvars.pivot(columns='LPERMNO', values='dlc')
dlcvars = dlcvars.fillna(value=None, method='ffill', limit=12)
txpvars = annualcompvars.pivot(columns='LPERMNO', values='txp')
txpvars = txpvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)

cf[cfmask] = ( (ibvars -
    ( ( (actvars - actvars.shift(12)) -
        (chevars - chevars.shift(12)) ) -
      ( (lctvars - lctvars.shift(12)) -
        (dlcvars - dlcvars.shift(12)) -
        (txpvars - txpvars.shift(12)) -
        (dpvars) ) ) ) /
      ( (atvars + atvars.shift(12)) / 2 ).replace(0,np.nan) )

cf.astype(float).to_csv(pathtoaddvars_1_200 + '3.csv')

#----------------------------------------------------------------------------------------------------------------------------------
	# 38 hire ( Employee Growth Rate )
'''hire=(emp-lag(emp))/lag(emp);
if missing(emp) or missing(lag(emp)) then hire=0
'''
empvars = annualcompvars.pivot(columns='LPERMNO', values='emp')
empvars = empvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

hire = (empvars - empvars.shift(12)) / empvars.shift(12)
hire[(empvars.isna()) & (empvars.shift(12)).isna()] = 0
hire = hire.replace(np.inf, 0)

hire = hire.replace(-np.inf, 0)

countvars[countvars.notnull()] = 1
hire *= countvars

hire.astype(float).to_csv(pathto_1_200 + '38.csv')

#----------------------------------------------------------------------------------------------------------------------------------
	# 84 ( Sales Growth )
'''sgr=(sale/lag(sale))-1
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)

sgr = ( salevars / salevars.shift(12).replace(0,np.nan) ) - 1

sgr.astype(float).to_csv(pathto_1_200 + '84.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 4 chpm
'''chpm = (ib / sale) - (lag(ib) / lag(sale))
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)

chpm = ( ibvars / salevars.replace(0,np.nan) ) - ( ibvars.shift(12) / salevars.shift(12).replace(0,np.nan) )

chpm.astype(float).to_csv(pathtoaddvars_1_200 + '4.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 5 chato
'''chato=(sale/((at+lag(at))/2)) - (lag(sale)/((lag(at)+lag2(at))/2));
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

chato=( salevars / ((atvars + atvars.shift(12)) / 2).replace(0,np.nan) ) - ( salevars.shift(12) / ((atvars.shift(12) + atvars.shift(24)) / 2).replace(0,np.nan) )

chato.astype(float).to_csv(pathtoaddvars_1_200 + '5.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 61 pchsale_pchinvt ( % Change in Sales - % Change in Inventory )
'''pchsale_pchinvt=((sale-lag(sale))/lag(sale))-((invt-lag(invt))/lag(invt))
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

pchsale_pchinvt=( (salevars - salevars.shift(12) ) / salevars.shift(12).replace(0,np.nan) ) - ( (invtvars - invtvars.shift(12) ) / invtvars.shift(12).replace(0,np.nan) )

pchsale_pchinvt.astype(float).to_csv(pathto_1_200 + '61.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 62 pchsale_pchrect ( % Change in Sales - % Change in A/R )
'''pchsale_pchrect=((sale-lag(sale))/lag(sale))-((rect-lag(rect))/lag(rect))
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
rectvars = annualcompvars.pivot(columns='LPERMNO', values='rect')
rectvars = rectvars.fillna(value=None, method='ffill', limit=12)

pchsale_pchrect = ( (salevars - salevars.shift(12) ) / salevars.shift(12).replace(0,np.nan) ) - ( (rectvars - rectvars.shift(12) ) / rectvars.shift(12).replace(0,np.nan) )

pchsale_pchrect.astype(float).to_csv(pathto_1_200 + '62.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 59 pchgm_pchsale ( % Change in Gross Margin - % Change in Sales )
'''pchgm_pchsale=(((sale-cogs)-(lag(sale)-lag(cogs)))/(lag(sale)-lag(cogs)))-((sale-lag(sale))/lag(sale))
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
cogsvars = annualcompvars.pivot(columns='LPERMNO', values='cogs')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=12)

pchgm_pchsale = ( ( ( salevars - cogsvars ) - ( salevars.shift(12) - cogsvars.shift(12) ) ) / ( salevars.shift(12) - cogsvars.shift(12) ).replace(0,np.nan) ) - ( ( salevars - salevars.shift(12) ) / salevars.shift(12).replace(0,np.nan) )

pchgm_pchsale.astype(float).to_csv(pathto_1_200 + '59.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 63 pchsale_pchxsga ( % Change in Sales - % Change in SG&A )
'''pchsale_pchxsga=( (sale-lag(sale))/lag(sale) )-( (xsga-lag(xsga)) /lag(xsga) )
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
xsgavars = annualcompvars.pivot(columns='LPERMNO', values='xsga')
xsgavars = xsgavars.fillna(value=None, method='ffill', limit=12)

pchsale_pchxsga = ( (salevars - salevars.shift(12) ) / salevars.shift(12).replace(0,np.nan) ) - ( (xsgavars - xsgavars.shift(12) ) / xsgavars.shift(12).replace(0,np.nan) )

pchsale_pchxsga.astype(float).to_csv(pathto_1_200 + '63.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 26 depr ( Depreciation / PP&E )
'''depr=dp/ppent
'''
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)

depr = dpvars / ppentvars.replace(0,np.nan)

depr.astype(float).to_csv(pathto_1_200 + '26.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 58 pchdepr ( % Change in Depreciation )
'''pchdepr=((dp/ppent)-(lag(dp)/lag(ppent)))/(lag(dp)/lag(ppent))
'''
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)

pchdepr = ( ( dpvars / ppentvars.replace(0,np.nan) ) - ( dpvars.shift(12) / ppentvars.shift(12).replace(0,np.nan) ) ) / ( dpvars.shift(12) / ppentvars.shift(12).replace(0,np.nan) ).replace(0,np.nan)

pchdepr.astype(float).to_csv(pathto_1_200 + '58.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 6 chadv
'''chadv=log(1+xad)-log((1+lag(xad)));	*had error here before, might work better now...
'''
xadvars = annualcompvars.pivot(columns='LPERMNO', values='xad')
xadvars = xadvars.fillna(value=None, method='ffill', limit=12)

chadv = np.log( xadvars + 1 ) - np.log( xadvars.shift(12) + 1 )

chadv.astype(float).to_csv(pathtoaddvars_1_200 + '6.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 42 invest ( Capital Expenditures and Inventory )
'''invest=( 	(ppegt-lag(ppegt)) +  (invt-lag(invt))	)	/ lag(at)
if missing(ppegt) then invest=( 	(ppent-lag(ppent)) +  (invt-lag(invt))	)	/ lag(at)
'''
ppegtvars = annualcompvars.pivot(columns='LPERMNO', values='ppegt')
ppegtvars = ppegtvars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
investmask = (ppegtvars.isna())

        # Calculating First Variable
invest = ( ( ppegtvars - ppegtvars.shift(12) ) + ( invtvars - invtvars.shift(12) ) ) / atvars.shift(12).replace(0,np.nan)

        # Applying Mask
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)

invest[investmask] = ( ( ppentvars - ppentvars.shift(12) ) + ( invtvars - invtvars.shift(12) ) ) / atvars.shift(12).replace(0,np.nan)

invest.astype(float).to_csv(pathto_1_200 + '42.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # 32 egr ( Growth in Common Shareholder Equity
'''egr=( (ceq-lag(ceq))/lag(ceq)  )
'''
ceqvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=12)

egr = ( ceqvars - ceqvars.shift(12) ) / ceqvars.shift(12).replace(0,np.nan)

egr.astype(float).to_csv(pathto_1_200 + '32.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variables 7 capx
'''if missing(capx) and count>=2 then capx=ppent-lag(ppent)
Use count variable
'''
capxvars = annualcompvars.pivot(columns='LPERMNO', values='capx')
capxvars = capxvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
capxmask = ( (capxvars.isna()) & (countvars >= 2) )

        # Applying Mask
capxvars[capxmask] = ppentvars - ppentvars.shift(12)

capx = capxvars.copy()

capx.astype(float).to_csv(pathtoaddvars_1_200 + '7.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 8 pchcapx
'''pchcapx = (capx - lag(capx)) / lag(capx)
Use additional Variable 7 ( capx )
'''
pchcapx = ( capx - capx.shift(12) ) / capx.shift(12).replace(0,np.nan)

pchcapx.astype(float).to_csv(pathtoaddvars_1_200 + '8.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 35 grCAPX ( Growth in Capital Expenditures )
'''grcapx = (capx - lag2(capx)) / lag2(capx)
Use additional Variable 7 ( capx )
'''
grcapx = ( capx - capx.shift(24) ) / capx.shift(24).replace(0,np.nan)

grcapx.astype(float).to_csv(pathto_1_200 + '35.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 9 grGW
'''grGW = (gdwl - lag(gdwl)) / lag(gdwl);
if missing(gdwl) or gdwl=0 then grGW=0;
if gdwl ne 0 and not missing(gdwl) and missing(grGW) then grGW=1;
'''
gdwlvars = annualcompvars.pivot(columns='LPERMNO', values='gdwl')
gdwlvars = gdwlvars.fillna(value=None, method='ffill', limit=12)

        # Calculating First Variable
grGW = ( gdwlvars - gdwlvars.shift(12) ) / gdwlvars.shift(12).replace(0,np.nan)

        # Creating Mask
grgwmask1 = ( (gdwlvars.isna()) | (gdwlvars==0) )
grgwmask2 = ( (gdwlvars!=0) & (gdwlvars.notnull()) & (grGW.isna()) )

        # Applying Mask
grGW[grgwmask1] = 0
grGW[grgwmask2] = 1

grGW.astype(float).to_csv(pathtoaddvars_1_200 + '9.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 10 woGW
'''if (not missing(gdwlia) and gdwlia ne 0) or (not missing(gdwlip) and gdwlip ne 0) or (not missing(gwo) and gwo ne 0) then woGW=1;
else woGW=0;
'''
gdwliavars = annualcompvars.pivot(columns='LPERMNO', values='gdwlia')
gdwliavars = gdwliavars.fillna(value=None, method='ffill', limit=12)
gdwlipvars = annualcompvars.pivot(columns='LPERMNO', values='gdwlip')
gdwlipvars = gdwlipvars.fillna(value=None, method='ffill', limit=12)
gwovars = annualcompvars.pivot(columns='LPERMNO', values='gwo')
gwovars = gwovars.fillna(value=None, method='ffill', limit=12)

        # Calculating First Variable
woGW = pd.DataFrame(columns=gdwliavars.columns, index=gdwliavars.index, data=0)

        # Creating Mask
wogwmask = ( ( (gdwliavars.notnull()) & (gdwliavars!=0) ) | ( (gdwlipvars.notnull()) & (gdwlipvars!=0) ) | ( (gwovars.notnull()) & (gwovars!=0) ) )

        # Applying Mask
woGW[wogwmask] = 1

woGW.astype(float).to_csv(pathtoaddvars_1_200 + '10.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 91 tang ( Dept Capacity / Firm Tangibility )
'''tang = (che + rect * 0.715 + invt * 0.547 + ppent * 0.535) / at
'''
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
rectvars = annualcompvars.pivot(columns='LPERMNO', values='rect')
rectvars = rectvars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

tang = (chevars + rectvars * 0.715 + invtvars * 0.547 + ppentvars * 0.535) / atvars.replace(0,np.nan)

tang.astype(float).to_csv(pathto_1_200 + '91.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 85 sin ( Sin Stocks )
'''if (2100 <= sic <= 2199) or (2080 <= sic <= 2085) or (
naics in ('7132', '71312', '713210', '71329', '713290', '72112', '721120'))
then
sin = 1; else sin = 0;
'''
sicvars = annualcompvars.pivot(columns='LPERMNO', values='sic')
sicvars = sicvars.fillna(value=None, method='ffill', limit=12)
naicsvars = annualcompvars.pivot(columns='LPERMNO', values='naics')
naicsvars = naicsvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
sinmask = ( ( ( 2100 <= sicvars) & (sicvars <= 2199) ) | ( (2080 <= sicvars) & (sicvars <= 2085) ) | ( naicsvars.isin(['7132', '71312', '713210', '71329', '713290', '72112', '721120']) ) )

        # Calculating First Variable
sin = pd.DataFrame(columns=sicvars.columns, index=sicvars.index, data=0)

        # Applying Mask
sin[sinmask] = 1

countvars[countvars.notnull()] = 1
sin *= countvars

sin.astype(float).to_csv(pathto_1_200 + '85.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 11 act
'''if missing(act) then act=che+rect+invt
'''
actvars = annualcompvars.pivot(columns='LPERMNO', values='act')
actvars = actvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
actmask = ( actvars.isna() )

        # Calculating First Variable
act = actvars.copy()

        # Applying Mask
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)
rectvars = annualcompvars.pivot(columns='LPERMNO', values='rect')
rectvars = rectvars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

act[actmask] = chevars + rectvars + invtvars

act.astype(float).to_csv(pathtoaddvars_1_200 + '11.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 12 lct
'''if missing(lct) then lct=ap
'''
lctvars = annualcompvars.pivot(columns='LPERMNO', values='lct')
lctvars = lctvars.fillna(value=None, method='ffill', limit=12)
apvars = annualcompvars.pivot(columns='LPERMNO', values='ap')
apvars = apvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
lctmask = ( lctvars.isna() )

        # Calculating First Variable
lct = lctvars.copy()

        # Applying Mask
lct[lctmask] = apvars

lct.astype(float).to_csv(pathtoaddvars_1_200 + '12.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 25 currat ( Current Ratio )
'''currat = act / lct
Use Additional Variable 11 ( act ) 879
Use Additional Variable 12 ( lct ) 902
'''
currat = act / lct.replace(0,np.nan)

currat.astype(float).to_csv(pathto_1_200 + '25.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 57 pchcurrat ( % Change in Current Ratio )
'''pchcurrat = ((act / lct) - (lag(act) / lag(lct))) / (lag(act) / lag(lct))
Use Additional Variable 11 ( act ) 879
Use Additional Variable 12 ( lct ) 902
'''
pchcurrat = ( ( act / lct.replace(0,np.nan) ) - ( act.shift(12) / lct.shift(12).replace(0,np.nan) ) ) / ( act.shift(12) / lct.shift(12).replace(0,np.nan) ).replace(0,np.nan)

pchcurrat.astype(float).to_csv(pathto_1_200 + '57.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 68 quick ( Quick Ratio )
'''quick = (act - invt) / lct
Use Additional Variable 11 ( act ) 879 
Use Additional Variable 12 ( lct ) 902
'''
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

quick = ( act - invtvars.astype(float) ) / lct.replace(0,np.nan)

quick.astype(float).to_csv(pathto_1_200 + '68.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 60 pchquick ( % Change in quick Ratio )
'''pchquick = (	(act -invt ) /lct - (lag(act ) -lag(invt) ) /lag(lct) )/  (   (lag(act) - lag(invt)) / lag(lct))
Use Additional Variable 9 ( act ) 879
Use Additional Variable 10 ( lct ) 902
'''
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

pchquick = (	(act - invtvars ) / lct.replace(0,np.nan) - ( act.shift(12)  - invtvars.shift(12) ) / lct.shift(12).replace(0,np.nan) ) /  (  ( act.shift(12) -  invtvars.shift(12) ) / lct.shift(12) ).replace(0,np.nan)

pchquick.astype(float).to_csv(pathto_1_200 + '60.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 79 salecash ( Sales to Cash )
'''salecash = sale / che
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
chevars = annualcompvars.pivot(columns='LPERMNO', values='che')
chevars = chevars.fillna(value=None, method='ffill', limit=12)

salecash = salevars / chevars.replace(0,np.nan)

salecash.astype(float).to_csv(pathto_1_200 + '79.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 81 salerec ( Sales to Receivables )
'''salerec = sale / rect
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
rectvars = annualcompvars.pivot(columns='LPERMNO', values='rect')
rectvars = rectvars.fillna(value=None, method='ffill', limit=12)

salerec = salevars / rectvars.replace(0,np.nan)

salerec.astype(float).to_csv(pathto_1_200 + '81.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 80 saleinv ( Sales to Inventory )
'''saleinv = sale / invt
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

saleinv = salevars / invtvars.replace(0,np.nan)

saleinv.astype(float).to_csv(pathto_1_200 + '80.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 64 pchsaleinv ( % Change Sales-to-Inventory )
'''pchsaleinv = ((sale / invt) - (lag(sale) / lag(invt))) / (lag(sale) / lag(invt))
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)

pchsaleinv = ( ( salevars / invtvars.replace(0,np.nan) ) - ( salevars.shift(12) / invtvars.shift(12).replace(0,np.nan) ) ) / ( salevars.shift(12) / invtvars.shift(12).replace(0,np.nan) ).replace(0,np.nan)

pchsaleinv.astype(float).to_csv(pathto_1_200 + '64.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 12 cashdebt ( Cash Flow to Debt )
'''cashdebt = (ib + dp) / ((lt + lag(lt)) / 2)
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)
ltvars = annualcompvars.pivot(columns='LPERMNO', values='lt')
ltvars = ltvars.fillna(value=None, method='ffill', limit=12)

cashdebt =  ( ibvars + dpvars ) / ( ( ltvars + ltvars.shift(12) )  / 2 ).replace(0,np.nan)

cashdebt.astype(float).to_csv(pathto_1_200 + '12.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 72 realestate ( Real Estate Holdings )
'''realestate = (fatb + fatl) / ppegt
if missing(ppegt) then realestate=(fatb+fatl) / ppent
'''
fatbvars = annualcompvars.pivot(columns='LPERMNO', values='fatb')
fatbvars = fatbvars.fillna(value=None, method='ffill', limit=12)
fatlvars = annualcompvars.pivot(columns='LPERMNO', values='fatl')
fatlvars = fatlvars.fillna(value=None, method='ffill', limit=12)
ppegtvars = annualcompvars.pivot(columns='LPERMNO', values='ppegt')
ppegtvars = ppegtvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)

fatbvars[fatbvars[fatlvars==0].isna()] = 0

        # Creating Mask
realestatemask = (ppegtvars.isna())

        # Calculating First Variable
realestate = ( fatbvars + fatlvars) / ppegtvars.replace(0,np.nan)

        # Applying Mask
realestate[realestatemask] = ( fatbvars + fatlvars) / ppentvars.replace(0,np.nan)

realestate.astype(float).to_csv(pathto_1_200 + '72.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 27 divi ( Dividend Initiation )
'''if (not missing(dvt) and dvt > 0) and (lag(dvt)=0 or missing(lag(dvt))) then divi=1; else divi=0
'''
dvtvars = annualcompvars.pivot(columns='LPERMNO', values='dvt')
dvtvars = dvtvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
divimask = ( ( (dvtvars.notnull()) & (dvtvars>0) ) & (dvtvars.shift(12)==0 ) | ( dvtvars.shift(12).isna() ) )

        # Calculating First Variable
divi = pd.DataFrame(columns=dvtvars.columns, index=dvtvars.index, data=0)

        # Applying Mask
divi[divimask] = 1

countvars[countvars.notnull()] = 1
divi *= countvars

divi.astype(float).to_csv(pathto_1_200 + '27.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 28 divo ( Dividen Omission )
'''if (missing(dvt) or dvt=0) and (lag(dvt) > 0 and not missing(lag(dvt))) then divo=1; else divo=0
'''
dvtvars = annualcompvars.pivot(columns='LPERMNO', values='dvt')
dvtvars = dvtvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
divomask = ( ( (dvtvars.isna()) | (dvtvars==0) ) & (dvtvars.shift(12)>0 ) & ( dvtvars.shift(12).notnull() ) )

        # Calculating First Variable
divo = pd.DataFrame(columns=dvtvars.columns, index=dvtvars.index, data=0)

        # Applying Mask
divo[divomask] = 1

countvars[countvars.notnull()] = 1
divo *= countvars

divo.astype(float).to_csv(pathto_1_200 + '28.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 13 obklg
'''obklg = ob / ((at + lag(at)) / 2)
'''
obvars = annualcompvars.pivot(columns='LPERMNO', values='ob')
obvars = obvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

obklg = obvars / ( atvars + atvars.shift(12) / 2 ).replace(0,np.nan)

obklg.astype(float).to_csv(pathtoaddvars_1_200 + '13.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 14 chobklg
'''chobklg = (ob - lag(ob)) / ((at + lag(at)) / 2)
'''
obvars = annualcompvars.pivot(columns='LPERMNO', values='ob')
obvars = obvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

chobklg = ( obvars - obvars.shift(12) ) / ( atvars + atvars.shift(12) / 2 ).replace(0,np.nan)

chobklg.astype(float).to_csv(pathtoaddvars_1_200 + '14.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 83 secureind ( Secured Debt Indicator )
'''if not missing(dm) and dm ne 0 then securedind=1; else securedind=0
'''
        # Creating Mask
secureindmask = ( (annualcompvars['dm'].notnull()) & (annualcompvars['dm']!=0) )

        # Calculating First Variable
annualcompvars['secureind'] = 0

        # Applying Mask
annualcompvars['secureind'][secureindmask] = 1

dmvars = annualcompvars.pivot(columns='LPERMNO', values='secureind')
dmvars = dmvars.fillna(value=None, method='ffill', limit=12)

securedind = dmvars.copy()

securedind.astype(float).to_csv(pathto_1_200 + '83.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 82 secured ( Secured Debt )
'''secured = dm / dltt
'''
dmvars = annualcompvars.pivot(columns='LPERMNO', values='dm')
dmvars = dmvars.fillna(value=None, method='ffill', limit=12)
dlttvars = annualcompvars.pivot(columns='LPERMNO', values='dltt')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=12)

secured = dmvars / dlttvars.replace(0,np.nan)

secured.astype(float).to_csv(pathto_1_200 + '82.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 24 convind ( convertible Debt Indicator )
'''if not missing(dc) and dc ne 0 or ( not missing(cshrc) and CSHRC ne 0) then convind=1; else convind=0
'''
    # Creating Mask
convindmask = ( ( (annualcompvars['dc'].notnull()) & (annualcompvars['dc']!=0) ) | ( (annualcompvars['cshrc'].notnull()) & (annualcompvars['cshrc']!=0) ) )

    # Calculating First Variable
annualcompvars['convind'] = 0

    # Applying Mask
annualcompvars['convind'][convindmask] = 1

convind = annualcompvars.pivot(columns='LPERMNO', values='convind')
convind = convind.fillna(value=None, method='ffill', limit=12)

convind.astype(float).to_csv(pathto_1_200 + '24.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 15 conv
'''conv = dc / dltt
'''
dcvars = annualcompvars.pivot(columns='LPERMNO', values='dc')
dcvars = dcvars.fillna(value=None, method='ffill', limit=12)
dlttvars = annualcompvars.pivot(columns='LPERMNO', values='dltt')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=12)

conv = dcvars / dlttvars.replace(0,np.nan)

conv.astype(float).to_csv(pathtoaddvars_1_200 + '15.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 36 grltnoa ( Growth in Long Term Net Operating Assets )
'''grltnoa = ((rect + invt + ppent + aco + intan + ao - ap - lco - lo) - 
(lag(rect) + lag(invt) + lag(ppent) + lag(aco) + lag(intan) + lag(ao) - lag(ap) - lag(lco) - lag(lo))- 
(rect - lag(rect) + invt - lag(invt) + aco - lag(aco) - (ap - lag(ap) + lco - lag(lco)) - dp)) / 
((at + lag(at)) / 2)
'''
rectvars = annualcompvars.pivot(columns='LPERMNO', values='rect')
rectvars = rectvars.fillna(value=None, method='ffill', limit=12)
invtvars = annualcompvars.pivot(columns='LPERMNO', values='invt')
invtvars = invtvars.fillna(value=None, method='ffill', limit=12)
ppentvars = annualcompvars.pivot(columns='LPERMNO', values='ppent')
ppentvars = ppentvars.fillna(value=None, method='ffill', limit=12)
acovars = annualcompvars.pivot(columns='LPERMNO', values='aco')
acovars = acovars.fillna(value=None, method='ffill', limit=12)
intanvars = annualcompvars.pivot(columns='LPERMNO', values='intan')
intanvars = intanvars.fillna(value=None, method='ffill', limit=12)
aovars = annualcompvars.pivot(columns='LPERMNO', values='ao')
aovars = aovars.fillna(value=None, method='ffill', limit=12)
apvars = annualcompvars.pivot(columns='LPERMNO', values='ap')
apvars = apvars.fillna(value=None, method='ffill', limit=12)
lcovars = annualcompvars.pivot(columns='LPERMNO', values='lco')
lcovars = lcovars.fillna(value=None, method='ffill', limit=12)
lovars = annualcompvars.pivot(columns='LPERMNO', values='lo')
lovars = lovars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)

grltnoa = ( ( ( (rectvars + invtvars + ppentvars + acovars + intanvars + aovars - apvars - lcovars - lovars) -
( rectvars.shift(12) + invtvars.shift(12) + ppentvars.shift(12) + acovars.shift(12) + intanvars.shift(12) + aovars.shift(12) - apvars.shift(12) - lcovars.shift(12) - lovars.shift(12) ) ) -
( rectvars - rectvars.shift(12) + invtvars - invtvars.shift(12) + acovars - acovars.shift(12) - ( apvars - apvars.shift(12) + lcovars - lcovars.shift(12) ) - dpvars ) ) /
( (atvars + atvars.shift(12) ) / 2).replace(0,np.nan) )

grltnoa.astype(float).to_csv(pathto_1_200 + '36.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 16 chdrc
'''chdrc = (dr - lag(dr)) / ((at + lag(at)) / 2)
'''
drvars = annualcompvars.pivot(columns='LPERMNO', values='dr')
drvars = drvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

chdrc = ( drvars - drvars.shift(12) ) / ( ( atvars + atvars.shift(12) ) / 2).replace(0,np.nan)

chdrc.astype(float).to_csv(pathtoaddvars_1_200 + '16.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 69 rd ( R&D Increase )
'''if ((xrd / at) - (lag(xrd / lag(at)))) / (lag(xrd / lag(at))) > .05 then rd=1; else rd=0
'''
xrdvars = annualcompvars.pivot(columns='LPERMNO', values='xrd')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
rdmask = ( (( (xrdvars / atvars.replace(0,np.nan)) - ( ( xrdvars / atvars.shift(12).replace(0,np.nan) ) ) ) / ( xrdvars / atvars.shift(12).replace(0,np.nan) )) > 0.5 )

        # Calculating First Variable
rd = pd.DataFrame(columns=xrdvars.columns, index=xrdvars.index, data=0)

        # Applying Mask
rd[rdmask] = 1

countvars[countvars.notnull()] = 1
rd *= countvars

rd.astype(float).to_csv(pathto_1_200 + '69.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variables 17 rdbias
'''rdbias = (xrd / lag(xrd)) - 1 - ib / lag(ceq)
'''
xrdvars = annualcompvars.pivot(columns='LPERMNO', values='xrd')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=12)
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
ceqvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=12)

rdbias = (xrdvars / xrdvars.shift(12).replace(0,np.nan)) - 1 - (ibvars / ceqvars.shift(12).replace(0,np.nan))

rdbias.astype(float).to_csv(pathtoaddvars_1_200 + '17.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 18 roe
'''roe = ib / lag(ceq)
'''
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
ceqvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=12)

roe = ibvars / ceqvars.shift(12).replace(0,np.nan)

roe.astype(float).to_csv(pathtoaddvars_1_200 + '18.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 54 operprof ( Operating Profitability )
'''operprof = (revt-cogs-xsga0-xint0)/lag(ceq)
'''
revtvars = annualcompvars.pivot(columns='LPERMNO', values='revt')
revtvars = revtvars.fillna(value=None, method='ffill', limit=12)
cogsvars = annualcompvars.pivot(columns='LPERMNO', values='cogs')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=12)
xsga0vars = annualcompvars.pivot(columns='LPERMNO', values='xsga0')
xsga0vars = xsga0vars.fillna(value=None, method='ffill', limit=12)
xint0vars = annualcompvars.pivot(columns='LPERMNO', values='xint0')
xint0vars = xint0vars.fillna(value=None, method='ffill', limit=12)
ceqvars = annualcompvars.pivot(columns='LPERMNO', values='ceq')
ceqvars = ceqvars.fillna(value=None, method='ffill', limit=12)

operprof = (revtvars-cogsvars-xsga0vars-xint0vars) / ceqvars.shift(12).replace(0,np.nan)
operprof.astype(float).to_csv(pathto_1_200 + '54.csv')
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
nivars = annualcompvars.pivot(columns='LPERMNO', values='ni')
nivars = nivars.fillna(value=None, method='ffill', limit=12)
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)
dlttvars = annualcompvars.pivot(columns='LPERMNO', values='dltt')
dlttvars = dlttvars.fillna(value=None, method='ffill', limit=12)
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)
cogsvars = annualcompvars.pivot(columns='LPERMNO', values='cogs')
cogsvars = cogsvars.fillna(value=None, method='ffill', limit=12)
scstkcvars = annualcompvars.pivot(columns='LPERMNO', values='scstkc')
scstkcvars = scstkcvars.fillna(value=None, method='ffill', limit=12)

ps = ( ( nivars [ nivars > 0 ] ) +
       ( oancfvars [ oancfvars > 0 ] ) +
       ( ( nivars / atvars.replace(0,np.nan) ) [ ( nivars / atvars.replace(0,np.nan) ) > ( nivars.shift(12) / atvars.shift(12).replace(0,np.nan) ) ] ) +
       ( oancfvars [ oancfvars > nivars ] ) +
       ( ( dlttvars / atvars.replace(0,np.nan) ) [ ( dlttvars / atvars.replace(0,np.nan) ) < ( dlttvars.shift(12) / atvars.shift(12).replace(0,np.nan) ) ] ) +
       ( ( actvars / lctvars.replace(0,np.nan) ) [ ( actvars / lctvars.replace(0,np.nan) ) > ( actvars.shift(12) / lctvars.shift(12).replace(0,np.nan) ) ] ) +
       ( ( ( salevars - cogsvars ) / salevars ) [ ( ( salevars - cogsvars ) / salevars ) > ( ( salevars.shift(12) - cogsvars.shift(12) ) / salevars.shift(12).replace(0,np.nan) )  ] ) +
       ( ( salevars / atvars.replace(0,np.nan) ) [ ( salevars / atvars )  > ( salevars.shift(12) / atvars.shift(12) ) ] ) +
       ( scstkcvars [ scstkcvars == 0 ] ) )

ps.astype(float).to_csv(pathto_1_200 + '67.csv')
#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 19 tb_1
''''*-----Lev and Nissim(2004);
if fyear <= 1978 then tr=.48;
if 1979 <= fyear <= 1986 then tr=.46;
if fyear=1987 then tr=.4;
if 1988 <= fyear <= 1992 then tr=.34;
if 1993 <= fyear then tr=.35;
tb_1 = ((txfo + txfed) / tr) / ib;
if missing(txfo) or missing(txfed) then tb_1=((txt-txdi) / tr) / ib; * they rank within industries;
if (txfo + txfed > 0 or txt > txdi) and ib <= 0 then tb_1=1;
'''
txfovars = annualcompvars.pivot(columns='LPERMNO', values='txfo')
txfovars = txfovars.fillna(value=None, method='ffill', limit=12)
txfedvars = annualcompvars.pivot(columns='LPERMNO', values='txfed')
txfedvars = txfedvars.fillna(value=None, method='ffill', limit=12)
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
txtvars = annualcompvars.pivot(columns='LPERMNO', values='txt')
txtvars = txtvars.fillna(value=None, method='ffill', limit=12)
txdivars = annualcompvars.pivot(columns='LPERMNO', values='txdi')
txdivars = txdivars.fillna(value=None, method='ffill', limit=12)
fyearvars = annualcompvars.pivot(columns='LPERMNO', values='fyear')
fyearvars = fyearvars.fillna(value=None, method='ffill', limit=12)

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

tb_1.astype(float).to_csv(pathtoaddvars_1_200 + '19.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 20 roa
'''roa = ni / ((at + lag(at)) / 2
For Mohanram score ( 2005 )
'''
nivars = annualcompvars.pivot(columns='LPERMNO', values='ni')
nivars = nivars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

roa = nivars / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

roa.astype(float).to_csv(pathtoaddvars_1_200 + '20.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 21 cfroa
'''cfroa = oancf / ((at + lag(at)) / 2);
if missing(oancf) then cfroa=(ib+dp) / ((at+lag(at)) / 2);
For Mohanram score ( 2005 )
'''
oancfvars = annualcompvars.pivot(columns='LPERMNO', values='oancf')
oancfvars = oancfvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

        # Creating Mask
cfroamask = (oancfvars.isna())

        # Calculating First Variable
cfroa = oancfvars / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

        # Applying Mask
ibvars = annualcompvars.pivot(columns='LPERMNO', values='ib')
ibvars = ibvars.fillna(value=None, method='ffill', limit=12)
dpvars = annualcompvars.pivot(columns='LPERMNO', values='dp')
dpvars = dpvars.fillna(value=None, method='ffill', limit=12)

cfroa[cfroamask] = ( ibvars + dpvars ) / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

cfroa.astype(float).to_csv(pathtoaddvars_1_200 + '21.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 22 xrdint
'''xrdint = xrd / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )
'''
xrdvars = annualcompvars.pivot(columns='LPERMNO', values='xrd')
xrdvars = xrdvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

xrdint = xrdvars / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

xrdint.astype(float).to_csv(pathtoaddvars_1_200 + '22.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 23 capxint
'''capxint = capx / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )
'''
capxvars = annualcompvars.pivot(columns='LPERMNO', values='capx')
capxvars = capxvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

capxint = capxvars / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

capxint.astype(float).to_csv(pathtoaddvars_1_200 + '23.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 24 xadint
'''xadint = xad / ((at + lag(at)) / 2)
For Mohanram score ( 2005 )
'''
xadvars = annualcompvars.pivot(columns='LPERMNO', values='xad')
xadvars = xadvars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)

xadint = xadvars / (( atvars + atvars.shift(12) ) / 2 ).replace(0,np.nan)

xadint.astype(float).to_csv(pathtoaddvars_1_200 + '24.csv')

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

cpi = pd.read_csv(pathtorawdata + 'cpi.csv')
dr = pd.date_range(start= ('1/1/' + str(cpi['Year'][0])), freq='M', periods=len(cpi) )
cpi.index = dr
cpi.index = pd.to_datetime(cpi.index).to_period('M').asfreq('D')
cpi = cpi['Value']

xsga0vars = annualcompvars.pivot(columns='LPERMNO', values='xsga0')
xsga0vars = xsga0vars.fillna(value=None, method='ffill', limit=12)
atvars = annualcompvars.pivot(columns='LPERMNO', values='at')
atvars = atvars.fillna(value=None, method='ffill', limit=12)
countvars = annualcompvars.pivot(columns='LPERMNO', values='count')
countvars = countvars.fillna(value=None, method='ffill', limit=12)

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

orgcap.to_csv(pathto_1_200 + '55.csv')

#----------------------------------------------------------------------------------------------------------------------------------

''' I am Skipping a part which I noted at the Misc File 233 - 241 '''

#----------------------------------------------------------------------------------------------------------------------------------
    # Description & Function for Industry Adjustments
''' Industry Adjustments for Variables '''
'''						/*other preparation steps for annual variables: industry adjustments*/
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
chpm = pd.read_csv('./data/AdditionalVariables/4.csv')

chpm.set_index('datadate', inplace=True)
chpm.index = pd.to_datetime(chpm.index).to_period('D')
chpm.columns = chpm.columns.map(int)

chpmia = indadj(chpm)

chpmia.astype(float).to_csv('./data/Anomalies/21.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 16 chatoia ( Industry-Adjusted Change in Asset Turnover )
''' See Description for Industry Adjustments
Use additional Variable 5 ( chato ) 
'''
chato = pd.read_csv('./data/AdditionalVariables/5.csv')

chato.set_index('datadate', inplace=True)
chato.index = pd.to_datetime(chato.index).to_period('D')
chato.columns = chato.columns.map(int)

chatoia = indadj(chato)

chatoia.astype(float).to_csv('./data/Anomalies/16.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 25 indsale
''' See Description for Industry Adjustments
'''
salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)

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

indsale.astype(float).to_csv(pathto_1_200 + '25.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 18 chempia ( Industry-Adjusted Change in Employees )
''' See Description for Industry Adjustments
Use Variable 38 ( hire )  
'''
hire = pd.read_csv('./data/Anomalies/38.csv')

hire.set_index('datadate', inplace=True)
hire.index = pd.to_datetime(hire.index).to_period('D')
hire.columns = hire.columns.map(int)

chempia = indadj(hire)

chempia.astype(float).to_csv('./data/Anomalies/18.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 10 bm_ia ( Industry-Adjusted Book to Market )
''' See Description for Industry Adjustments
Use Variable 9 ( bm ) 
'''
bm = pd.read_csv('./data/Anomalies/9.csv')

bm.set_index('datadate', inplace=True)
bm.index = pd.to_datetime(bm.index).to_period('D')
bm.columns = bm.columns.map(int)


bm_ia = indadj(bm)

bm_ia.astype(float).to_csv('./data/Anomalies/10.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 56 pchcapx_ia ( Industry-Adjusted % chang in Capital Expenditures )
'''See Description for Industry Adjustments
Use Additional Variable 8 ( pchcapx ) 
'''
pchcapx = pd.read_csv('./data/AdditionalVariables/8.csv')

pchcapx.set_index('datadate', inplace=True)
pchcapx.index = pd.to_datetime(pchcapx.index).to_period('D')
pchcapx.columns = pchcapx.columns.map(int)

pchcapx_ia = indadj(pchcapx)

pchcapx_ia.astype(float).to_csv('./data/Anomalies/56.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 92 tb ( Tax Income to Book Income )
'''See Description for Industry Adjustments
Use Additional Variable 19 ( tb_1 ) 
'''
tb_1 = pd.read_csv('./data/AdditionalVariables/19.csv')

tb_1.set_index('datadate', inplace=True)
tb_1.index = pd.to_datetime(tb_1.index).to_period('D')
tb_1.columns = tb_1.columns.map(int)

tb = indadj(tb_1)

tb.astype(float).to_csv('./data/Anomalies/92.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 15 cfp_ia ( Industry_Adjusted Cash Flow to Price Ratio )
'''See Description for Industry Adjustments
Use Variable 14 ( cfp ) 
'''
cfp = pd.read_csv('./data/Anomalies/14.csv')

cfp.set_index('datadate', inplace=True)
cfp.index = pd.to_datetime(cfp.index).to_period('D')
cfp.columns = cfp.columns.map(int)

cfp_ia = indadj(cfp)

cfp_ia.astype(float).to_csv('./data/Anomalies/15.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 52 mve_ia ( Industry_Adjusted size )
'''See Description for Industry Adjustments
'''
mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

mve_ia = indadj(mve_f)

mve_ia.astype(float).to_csv('./data/Anomalies/52.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    #  37 herf ( Industry Sale Concentration  ) 
'''See Description for Industry Adjustments
sum( (sale/indsale)*(sale/indsale) ) as herf
'''
indsale = pd.read_csv(pathtoaddvars_1_200 + '25.csv')

indsale.set_index('datadate', inplace=True)
indsale.index = pd.to_datetime(indsale.index).to_period('D')
indsale.columns = indsale.columns.map(int)

salevars = annualcompvars.pivot(columns='LPERMNO', values='sale')
salevars = salevars.fillna(value=None, method='ffill', limit=12)

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

herf.astype(float).to_csv(pathto_1_200 + '37.csv')


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
m4[m4.isna() & m3.notnull()] = .5

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
m6[m6.isna() & m3.notnull()] = .5
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
m7[m7.isna() & m3.notnull()] = .5
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
m8[m8.isna() & m3.notnull()] = .5

ms = m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8

countvars = monthlycrspvars.pivot(columns='PERMNO', values='count')
countvars[countvars.notnull()] = 1

ms = ms*countvars
ms.astype(float).to_csv(pathto_1_200 + '50.csv')

#----------------------------------------------------------------------------------------------------------------------------------


# =============================================================================
# QUARTERLY VARIABLES
# =============================================================================


#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 26 pstk
'''if not missing(pstkrq) then pstk=pstkrq;
else pstk=pstkq
'''
pstkrqvars = quarterlycompvars.pivot(columns='LPERMNO', values='pstkrq')
pstkrqvars = pstkrqvars.fillna(value=None, method='ffill', limit=3)
pstkqvars = quarterlycompvars.pivot(columns='LPERMNO', values='pstkq')
pstkqvars = pstkqvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
pstkmask1 = (pstkrqvars.notnull())
pstkmask2 = (pstkrqvars.isna())

        # Calculating First Variable
pstk = pd.DataFrame(index=pstkrqvars.index, columns=pstkrqvars.columns)

        # Applying Mask
pstk[pstkmask1] = pstkrqvars.copy()
pstk[pstkmask2] = pstkqvars.copy()

pstk.astype(float).to_csv('./data/AdditionalVariables/26.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 27 scal
'''scal=seqq
if missing(seqq) then scal=ceqq+pstk;
if missing(seqq) and (missing(ceqq) or missing(pstk)) then scal=atq-ltq
Use Additional Variable 26 ( pstk )
'''
seqqvars = quarterlycompvars.pivot(columns='LPERMNO', values='seqq')
seqqvars = seqqvars.fillna(value=None, method='ffill', limit=3)
ceqqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ceqq')
ceqqvars = ceqqvars.fillna(value=None, method='ffill', limit=3)
atqvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atqvars = atqvars.fillna(value=None, method='ffill', limit=3)
ltqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ltq')
ltqvars = ltqvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
scalmask1 = (seqqvars.isna())
scalmask2 = ( (seqqvars.isna()) & ( (ceqqvars.isna()) | (pstk.isna()) ) )

        # Calculating First Variable
scal = seqqvars.copy()

        # Applying Mask
scal[scalmask1] = ceqqvars + pstk
scal[scalmask2] = atqvars - ltqvars

scal.astype(float).to_csv('./data/AdditionalVariables/27.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 22 chtx ( Change in Tax Expense )
'''chtx=(txtq-lag4(txtq))/lag4(atq)
'''
txtqvars = quarterlycompvars.pivot(columns='LPERMNO', values='txtq')
txtqvars = txtqvars.fillna(value=None, method='ffill', limit=3)
atqvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atqvars = atqvars.fillna(value=None, method='ffill', limit=3)

chtx=( txtqvars - txtqvars.shift(12) ) / atqvars.replace(0,np.nan).shift(12)

chtx.astype(float).to_csv('./data/Anomalies/22.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 74 roaq ( Return on Assets )
'''roaq=ibq/lag(atq)
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibqvars = ibqvars.fillna(value=None, method='ffill', limit=3)
atqvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atqvars = atqvars.fillna(value=None, method='ffill', limit=3)


roaq=ibqvars / atqvars.replace(0,np.nan).shift(3)

roaq.astype(float).to_csv('./data/Anomalies/74.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 76 roeq ( Return on Equity )
'''roeq=(ibq)/lag(scal)
Use Additional Variable 27 ( scal )
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibqvars = ibqvars.fillna(value=None, method='ffill', limit=3)

roeq = ibqvars / scal.shift(3).replace(0,np.nan)

roeq.astype(float).to_csv('./data/Anomalies/76.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 78 rsup ( Revenue Surprise )
'''rsup=(saleq-lag4(saleq))/mveq
'''
saleqvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
saleqvars = saleqvars.fillna(value=None, method='ffill', limit=3)

mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')

rsup = ( saleqvars - saleqvars.shift(12) ) / mve_f.replace(0,np.nan)

rsup.astype(float).to_csv('./data/Anomalies/78.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 28 sacc
'''sacc=( (actq-lag(actq) - (cheq-lag(cheq))) - (  (lctq-lag(lctq))-(dlcq-lag(dlcq)) ) ) /saleq; ;
  if saleq<=0 then sacc=( (actq-lag(actq) - (cheq-lag(cheq))) - (  (lctq-lag(lctq))-(dlcq-lag(dlcq)) ) ) /.01
'''
actqvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
actqvars = actqvars.fillna(value=None, method='ffill', limit=3)
cheqvars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
cheqvars = cheqvars.fillna(value=None, method='ffill', limit=3)
lctqvars = quarterlycompvars.pivot(columns='LPERMNO', values='lctq')
lctqvars = lctqvars.fillna(value=None, method='ffill', limit=3)
dlcqvars = quarterlycompvars.pivot(columns='LPERMNO', values='dlcq')
dlcqvars = dlcqvars.fillna(value=None, method='ffill', limit=3)
saleqvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
saleqvars = saleqvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
saccmask = ( saleqvars <= 0 )

        # Calculating First Variable
sacc = ( ( ( actqvars - actqvars.shift(3) ) -
         ( cheqvars - cheqvars.shift(3) ) ) -
       ( ( lctqvars - lctqvars.shift(3) ) -
         ( dlcqvars - dlcqvars.shift(3) ) ) ) / saleqvars

        # Applying Mask
sacc[saccmask] = ( ( ( actqvars - actqvars.shift(3) ) -
                     ( cheqvars - cheqvars.shift(3) ) ) -
                   ( ( lctqvars - lctqvars.shift(3) ) -
                     ( dlcqvars - dlcqvars.shift(3) ) ) ) / .01

sacc.astype(float).to_csv('./data/AdditionalVariables/28.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 89 stdacc ( Accrual Volatility )
'''stdacc = std( sacc,lag(sacc),lag2(sacc),lag3(sacc),lag4(sacc),lag5(sacc),lag6(sacc),lag7(sacc),
lag8(sacc),lag9(sacc),lag10(sacc),lag11(sacc),lag12(sacc),lag13(sacc),lag14(sacc),lag15(sacc))
Use Additional Variable 28 ( sacc )
'''
actqvars = quarterlycompvars.pivot(columns='LPERMNO', values='actq')
a = actqvars.copy()
a[a.notnull()] = 1

sacc2 = sacc * a

stdacc = sacc2.astype(float).rolling(window=45, min_periods=12).std()
stdacc = stdacc[sacc2.notnull()]

stdacc = stdacc.fillna(value=None, method='ffill', limit=3)

stdaccmask = pd.DataFrame(np.ones(shape=sacc.shape), index=sacc.index, columns=sacc.columns)
stdaccmask[sacc.isna()] = np.nan

stdacc *= stdaccmask

stdacc.astype(float).to_csv(pathto_1_200 + '89.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 29 sgrvol
'''sgrvol=std(rsup,lag(rsup),lag2(rsup),lag3(rsup),lag4(rsup),lag5(rsup),lag6(rsup),lag7(rsup),
lag8(rsup),lag9(rsup),lag10(rsup),lag11(rsup),lag12(rsup),lag13(rsup),lag14(rsup))
Use Variable 78 ( rsup )
'''

saleqvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
a = saleqvars.copy()
a[a.notnull()] = 1

rsup2 = rsup * a

sgrvol = rsup2.astype(float).rolling(window=42, min_periods=11).std()
sgrvol = sgrvol[rsup2.notnull()]

sgrvol = sgrvol.fillna(value=None, method='ffill', limit=3)

sgrvolmask = pd.DataFrame(np.ones(shape=rsup.shape), index=rsup.index, columns=rsup.columns)
sgrvolmask[rsup.isna()] = np.nan

sgrvol *= sgrvolmask

sgrvol.astype(float).to_csv(pathtoaddvars_1_200 + '29.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 75 roavol ( Earnings Volatility )
'''roavol=std(roaq,lag(roaq),lag2(roaq),lag3(roaq),lag4(roaq),lag5(roaq),lag6(roaq),lag7(roaq),
lag8(roaq),lag9(roaq),lag10(roaq),lag11(roaq),lag12(roaq),lag13(roaq),lag14(roaq),lag15(roaq))
Use Additional Variable 74 ( roaq )
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
a = ibqvars.copy()
a[a.notnull()] = 1

roaq2 = roaq* a

roavol = roaq2.astype(float).rolling(window=45, min_periods=12).std()
roavol = roavol[roaq2.notnull()]

roavol = roavol.fillna(value=None, method='ffill', limit=3)

roavolmask = pd.DataFrame(np.ones(shape=roaq.shape), index=roaq.index, columns=roaq.columns)
roavolmask[roaq.isna()] = np.nan

roavol *= roavolmask

roavol.astype(float).to_csv(pathto_1_200 + '75.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 30 ( scf )
'''scf=(ibq/saleq)-sacc;
if saleq<=0 then scf=(ibq/.01)-sacc
Use Additional Variable 27 sacc lines between 1634 - 1663
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibqvars = ibqvars.fillna(value=None, method='ffill', limit=3)
saleqvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
saleqvars = saleqvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
scfmask = (saleqvars <= 0)

        # Calculating First Variable
scf = ( ibqvars / saleqvars ) - sacc

        # Applying Mask
scf[scfmask] = ( ibqvars / .01 ) - sacc

scf.astype(float).to_csv('./data/AdditionalVariables/30.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 90 stdcf ( Cash Flow Volatility )
'''stdcf=std(scf,lag(scf),lag2(scf),lag3(scf),lag4(scf),lag5(scf),lag6(scf),lag7(scf),
lag8(scf),lag9(scf),lag10(scf),lag11(scf),lag12(scf),lag13(scf),lag14(scf),lag15(scf))
Use Additional Variable 29 scf lines between 1716 - 1735
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
a = ibqvars.copy()
a[a.notnull()] = 1

scf2 = scf * a

stdcf = scf2.astype(float).rolling(window=45, min_periods=15).std()
stdcf = stdcf[scf2.notnull()]

stdcf = stdcf.fillna(value=None, method='ffill', limit=3)

stdcf.astype(float).to_csv('./data/Anomalies/90.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 11 cash ( Cash Holdings )
'''cash=cheq/atq
'''
cheqvars = quarterlycompvars.pivot(columns='LPERMNO', values='cheq')
cheqvars = cheqvars.fillna(value=None, method='ffill', limit=3)
atqvars = quarterlycompvars.pivot(columns='LPERMNO', values='atq')
atqvars = atqvars.fillna(value=None, method='ffill', limit=3)


cash = cheqvars / atqvars.replace(0,np.nan)

cash.astype(float).to_csv('./data/Anomalies/11.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 23 cinvest ( Corporate Investments )
'''cinvest=((ppentq-lag(ppentq))/saleq)-mean(((lag(ppentq)-lag2(ppentq))/lag(saleq)),((lag2(ppentq)-lag3(ppentq))/lag2(saleq)),((lag3(ppentq)-lag4(ppentq))/lag3(saleq)));
if saleq<=0 then cinvest=((ppentq-lag(ppentq))/.01)-mean(((lag(ppentq)-lag2(ppentq))/(.01)),((lag2(ppentq)-lag3(ppentq))/(.01)),((lag3(ppentq)-lag4(ppentq))/(.01)))
'''
ppentqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ppentq')
ppentqvars = ppentqvars.fillna(value=None, method='ffill', limit=3)
saleqvars = quarterlycompvars.pivot(columns='LPERMNO', values='saleq')
saleqvars = saleqvars.fillna(value=None, method='ffill', limit=3)

        # Creating Mask
cinvestmask = (saleqvars<=0)

        # Calculating First Variable
cinvest=( ( ( ppentqvars - ppentqvars.shift(3) ) / saleqvars.replace(0,np.nan) )-
        ( ( ( ppentqvars.shift(3) - ppentqvars.shift(6) / saleqvars.replace(0,np.nan).shift(3) ) ) +
          ( ( ppentqvars.shift(6) - ppentqvars.shift(9) / saleqvars.replace(0,np.nan).shift(6) ) ) +
          ( ( ppentqvars.shift(9) - ppentqvars.shift(12) / saleqvars.replace(0,np.nan).shift(9) ) ) ) / 3 )

        # Appyling Mask
cinvest[cinvestmask] = ( ( ( ppentqvars - ppentqvars.shift(3) ) / .01 )-
                       ( ( ( ppentqvars.shift(3) - ppentqvars.shift(6) / .01 ) ) +
                         ( ( ppentqvars.shift(6) - ppentqvars.shift(9) / .01 ) ) +
                         ( ( ppentqvars.shift(9) - ppentqvars.shift(12) / .01 ) ) ) / 3 )

cinvest.astype(float).to_csv('./data/Anomalies/23.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 31 che
'''che=ibq-lag4(ibq)
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibqvars = ibqvars.fillna(value=None, method='ffill', limit=3)

che = ibqvars - ibqvars.shift(12)

che.astype(float).to_csv('./data/AdditionalVariables/31.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 53 nincr ( Number of Earnings Increases )
'''nincr	=(  (ibq>lag(ibq)) + 
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq)) +
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq)) + 
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq)) + 
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq)) + 
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq)) + 
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq))*(lag6(ibq)>lag7(ibq)) +
(ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq))*(lag6(ibq)>lag7(ibq))*(lag7(ibq)>lag8(ibq))  )
'''
ibqvars = quarterlycompvars.pivot(columns='LPERMNO', values='ibq')
ibqvars = ibqvars.fillna(value=None, method='ffill', limit=3)

nincr	= ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] )  * ( ( ibqvars.shift(9) ) [ ibqvars.shift(9) > ibqvars.shift(12) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] )  * ( ( ibqvars.shift(9) ) [ ibqvars.shift(9) > ibqvars.shift(12) ] )  * ( ( ibqvars.shift(12) ) [ ibqvars.shift(12) > ibqvars.shift(15) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] )  * ( ( ibqvars.shift(9) ) [ ibqvars.shift(9) > ibqvars.shift(12) ] )  * ( ( ibqvars.shift(12) ) [ ibqvars.shift(12) > ibqvars.shift(15) ] )  * ( ( ibqvars.shift(15) ) [ ibqvars.shift(15) > ibqvars.shift(18) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] )  * ( ( ibqvars.shift(9) ) [ ibqvars.shift(9) > ibqvars.shift(12) ] )  * ( ( ibqvars.shift(12) ) [ ibqvars.shift(12) > ibqvars.shift(15) ] )  * ( ( ibqvars.shift(15) ) [ ibqvars.shift(15) > ibqvars.shift(18) ] )  * ( ( ibqvars.shift(18) ) [ ibqvars.shift(18) > ibqvars.shift(21) ] ) ) +
          ( ( ibqvars [ ibqvars > ibqvars.shift(3) ] ) * ( ( ibqvars.shift(3) ) [ ibqvars.shift(3) > ibqvars.shift(6) ] )  * ( ( ibqvars.shift(6) ) [ ibqvars.shift(6) > ibqvars.shift(9) ] )  * ( ( ibqvars.shift(9) ) [ ibqvars.shift(9) > ibqvars.shift(12) ] )  * ( ( ibqvars.shift(12) ) [ ibqvars.shift(12) > ibqvars.shift(15) ] )  * ( ( ibqvars.shift(15) ) [ ibqvars.shift(15) > ibqvars.shift(18) ] )  * ( ( ibqvars.shift(18) ) [ ibqvars.shift(18) > ibqvars.shift(21) ] )  * ( ( ibqvars.shift(21) ) [ibqvars.shift(21) > ibqvars.shift(8) ] ) ) )

nincr.astype(float).to_csv('./data/Anomalies/53.csv')

#----------------------------------------------------------------------------------------------------------------------------------
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

aeavol = aeavol.fillna(value=None, method='ffill', limit=3)

aeavol.to_csv(pathto_1_200 + '3.csv')

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
ear = ear.fillna(value=None, method='ffill', limit=3)

ear.to_csv(pathto_1_200 + '31.csv')

#----------------------------------------------------------------------------------------------------------------------------------

'''There are some other variables that are based on monthly-CRSP information (already in the dataset from monthly CRSP)
create those variables plus a couple of others'''



#----------------------------------------------------------------------------------------------------------------------------------
    # 49 mom6m ( 6-Month Momentum )
'''mom6m=  (  (1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret)) ) - 1
'''
mom6mvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
mom6mvars[mom6mvars<0] = np.nan
mom6mvars = mom6mvars.pct_change(fill_method=None)

mom6m =  (  ( 1 + mom6mvars.shift(2) )  * ( 1 + mom6mvars.shift(3) ) * ( 1 + mom6mvars.shift(4) ) * ( 1 + mom6mvars.shift(5) ) * ( 1 + mom6mvars.shift(6) ) ) - 1

mom6m.astype(float).to_csv(pathto_1_200 + '49.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 46 mom12m ( 12-Month Momentum )
'''mom12m=  (   (1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret))*
(1+lag7(ret))*(1+lag8(ret))*(1+lag9(ret))*(1+lag10(ret))*(1+lag11(ret))*(1+lag12(ret))   ) - 1
'''
mom12mvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
mom12mvars[mom12mvars<0] = np.nan
mom12mvars = mom12mvars.pct_change(fill_method=None)

mom12m =  (  ( 1 + mom12mvars.shift(2) )  * ( 1 + mom12mvars.shift(3) ) * ( 1 + mom12mvars.shift(4) ) * ( 1 + mom12mvars.shift(5) ) * ( 1 + mom12mvars.shift(6) ) *
             ( 1 + mom12mvars.shift(7) )  * ( 1 + mom12mvars.shift(8) ) * ( 1 + mom12mvars.shift(9) ) * ( 1 + mom12mvars.shift(10) ) * ( 1 + mom12mvars.shift(11) ) *
             ( 1 + mom12mvars.shift(12) ) ) - 1

mom12m.astype(float).to_csv('./data/Anomalies/46.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 48 mom36m ( 36-Month Momentum )
'''mom36m=(   (1+lag13(ret))*(1+lag14(ret))*(1+lag15(ret))*(1+lag16(ret))*(1+lag17(ret))*(1+lag18(ret))   *
(1+lag19(ret))*(1+lag20(ret))*(1+lag21(ret))*(1+lag22(ret))*(1+lag23(ret))*(1+lag24(ret))*
(1+lag25(ret))*(1+lag26(ret))*(1+lag27(ret))*(1+lag28(ret))*(1+lag29(ret))*(1+lag30(ret))     *
(1+lag31(ret))*(1+lag32(ret))*(1+lag33(ret))*(1+lag34(ret))*(1+lag35(ret))*(1+lag36(ret))  ) - 1
'''
mom36mvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
mom36mvars[mom36mvars<0] = np.nan
mom36mvars = mom36mvars.pct_change(fill_method=None)

mom36m =  (  ( 1 + mom36mvars.shift(13) ) * ( 1 + mom36mvars.shift(14) ) * ( 1 + mom36mvars.shift(15) ) * ( 1 + mom36mvars.shift(16) ) * ( 1 + mom36mvars.shift(17) ) *
             ( 1 + mom36mvars.shift(18) ) * ( 1 + mom36mvars.shift(19) ) * ( 1 + mom36mvars.shift(20) ) * ( 1 + mom36mvars.shift(21) ) * ( 1 + mom36mvars.shift(22) ) *
             ( 1 + mom36mvars.shift(23) ) * ( 1 + mom36mvars.shift(24) ) * ( 1 + mom36mvars.shift(25) ) * ( 1 + mom36mvars.shift(26) ) * ( 1 + mom36mvars.shift(27) ) *
             ( 1 + mom36mvars.shift(28) ) * ( 1 + mom36mvars.shift(29) ) * ( 1 + mom36mvars.shift(30) ) * ( 1 + mom36mvars.shift(31) ) * ( 1 + mom36mvars.shift(32) ) *
             ( 1 + mom36mvars.shift(33) ) * ( 1 + mom36mvars.shift(34) ) * ( 1 + mom36mvars.shift(35) ) * ( 1 + mom36mvars.shift(36) ) ) - 1

mom36m.astype(float).to_csv('./data/Anomalies/48.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 47 mom1m ( 1-Month Momentum )
'''mom1m=	lag(ret)
'''
mom1mvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
mom1mvars[mom1mvars<0] = np.nan
mom1mvars = mom1mvars.pct_change(fill_method=None)

mom1m = mom1mvars.shift()

mom1m.astype(float).to_csv('./data/Anomalies/47.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 29 dolvol ( Dollar Trading Volume )
'''dolvol=log(lag2(vol)*lag2(prc))
'''
prcvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
prcvars[prcvars<0] = np.nan
volvars = monthlycrspvars.pivot(columns='PERMNO', values='VOL')

dolvol = np.log( (volvars.shift(2) * prcvars.shift(2)).replace(0,np.nan) )

dolvol.astype(float).to_csv('./data/Anomalies/29.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 20 chmom ( Change in 6-Month Momentum )
'''chmom =(   (1+lag(ret))*(1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret))   ) - 1
- ((  (1+lag7(ret))*(1+lag8(ret))*(1+lag9(ret))*(1+lag10(ret))*(1+lag11(ret))*(1+lag12(ret))   ) - 1)
'''
chmomvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
chmomvars[chmomvars<0] = np.nan
chmomvars = chmomvars.pct_change(fill_method=None)

chmom = ( ( ( ( 1 + chmomvars.shift() ) * ( 1 + chmomvars.shift(2) ) * ( 1 + chmomvars.shift(3) ) * ( 1 + chmomvars.shift(4) ) * ( 1 + chmomvars.shift(5) ) * ( 1 + chmomvars.shift(6) ) ) - 1 ) -
        ( ( ( 1 + chmomvars.shift(7) ) * ( 1 + chmomvars.shift(8) ) * ( 1 + chmomvars.shift(9) ) * ( 1 + chmomvars.shift(10) ) * ( 1 + chmomvars.shift(11) ) * ( 1 + chmomvars.shift(12) ) ) - 1) )
chmom.astype(float).to_csv('./data/Anomalies/20.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 93 turn ( Share Turnover )
'''turn=mean(lag(vol),lag2(vol),lag3(vol))/shrout
'''
volvars = monthlycrspvars.pivot(columns='PERMNO', values='VOL')
shroutvars = monthlycrspvars.pivot(columns='PERMNO', values='SHROUT')

vol1 = volvars.shift(1)
vol2 = volvars.shift(2)
vol3 = volvars.shift(3)

turn = ( ( vol1 + vol2 + vol3 ) / 3 ) / shroutvars.replace(0,np.nan)
turn.astype(float).to_csv('./data/Anomalies/93.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 31 retcons_pos
'''if lag(ret)>0 and lag2(ret)>0 and lag3(ret)>0 and lag4(ret)>0 and lag5(ret)>0 and lag6(ret)>0 then retcons_pos=1; else retcons_pos=0
'''
retcons_posvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
retcons_posvars[retcons_posvars<0] = np.nan
retcons_posvars = retcons_posvars.pct_change(fill_method=None)

        # Creating Mask
retcons_posmask = ( (retcons_posvars.shift() > 0) & (retcons_posvars.shift(2) > 0) & (retcons_posvars.shift(3) > 0) & (retcons_posvars.shift(4) > 0) & (retcons_posvars.shift(5) > 0) & (retcons_posvars.shift(6) > 0) )

        # Calculating First Variable
retcons_pos = pd.DataFrame(index=retcons_posvars.index, columns=retcons_posvars.columns, data=0)

        # Applying Mask
retcons_pos[retcons_posmask] = 1

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 32 retcons_neg
'''if lag(ret)<0 and lag2(ret)<0 and lag3(ret)<0 and lag4(ret)<0 and lag5(ret)<0 and lag6(ret)<0 then retcons_neg=1; else retcons_neg=0
'''
retcons_negvars = monthlycrspvars.pivot(columns='PERMNO', values='ALTPRC')
retcons_negvars[retcons_negvars<0] = np.nan
retcons_negvars = retcons_negvars.pct_change(fill_method=None)

        # Creating Mask
retcons_negmask = ( (retcons_negvars.shift() < 0) & (retcons_negvars.shift(2) < 0) & (retcons_negvars.shift(3) < 0) & (retcons_negvars.shift(4) < 0) & (retcons_negvars.shift(5) < 0) & (retcons_negvars.shift(6) < 0) )

        # Calculating First Variable
retcons_neg = pd.DataFrame(index=retcons_negvars.index, columns=retcons_negvars.columns, data=0)

        # Applying Mask
retcons_neg[retcons_negmask] = 1

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 33 IPO
'''if count<=12 then IPO=1; else IPO=0
'''
countvars = monthlycrspvars.pivot(columns='PERMNO', values='count')

        # Creating Mask
ipomask = ( countvars <= 12 )

        # Calculating First Variable
ipo = pd.DataFrame(index=countvars.index, columns=countvars.columns, data=0)

        # Applying Mask
ipo[ipomask] = 1

#----------------------------------------------------------------------------------------------------------------------------------
    # 41 indmom ( Industry Momentum )
'''mean(mom12m) as indmom
Use Variable 46 mom12m ( 12-Month Momentum )
'''
sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')

mom12m = pd.read_csv('./data/Anomalies/46.csv')
mom12m.set_index('date', inplace=True)
mom12m.index = pd.to_datetime(mom12m.index).to_period('D')
mom12m.columns = mom12m.columns.map(int)

indmom = pd.DataFrame(index=mom12m.index, columns=mom12m.columns)

for i in indmom.index:
    a = mom12m.loc[i].dropna()
    b = sic.loc[i]
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

    c.loc[d.index, 'values'] = d
    c.set_index('PERMNO', inplace=True)

    indmom.loc[i,c.index] = c['values']

indmom.astype(float).to_csv('./data/Anomalies/41.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' 501 - SIC Codes
'''
    # 501 sic ( SIC Codes )
sic = monthlycrspvars.pivot(columns='PERMNO', values='SICCD2')

sic = sic.replace('na', np.nan)
sic = sic.replace('Z', np.nan)

sic.astype(float).to_csv('./data/Anomalies500-600/501.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' Monthly Logarithmic Market Cap
'''
mve_f = monthlycrspvars.pivot(columns='PERMNO', values='mve_f')
mvel1 = np.log(mve_f)

mvel1.astype(float).to_csv(pathto_1_200 + '51.csv')

#----------------------------------------------------------------------------------------------------------------------------------
''' finally, a few more directly from daily CRSP to create monthly variables '''
#----------------------------------------------------------------------------------------------------------------------------------
    # 45 maxret ( Maximum Daily Return )
'''max(ret) as maxret
'''
retvars = pd.read_csv(pathtorawdata + 'dailyret.csv')
retvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
retvars = retvars.pivot(index='date', columns='PERMNO', values='RET')
retvars.index = pd.to_datetime(retvars.index, format='%d/%m/%Y')
retvars = retvars.sort_index()

retvars = retvars.replace('C', np.nan)
retvars = retvars.replace('B', np.nan)
retvars = retvars.astype('float')

maxretvars = retvars.copy()
maxretvars.index = pd.to_datetime(maxretvars.index).to_period('M')

maxret = maxretvars.groupby(by=maxretvars.index).max()
maxret.index = maxret.index.asfreq('D')

maxret.astype(float).to_csv(pathto_1_200 + '45.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 73 retvol ( Return Volatility )
'''std(ret) as retvol
'''
retvars = pd.read_csv(pathtorawdata + 'dailyret.csv')
retvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
retvars = retvars.pivot(index='date', columns='PERMNO', values='RET')
retvars.index = pd.to_datetime(retvars.index, format='%d/%m/%Y')
retvars = retvars.sort_index()

retvars = retvars.replace('C', np.nan)
retvars = retvars.replace('B', np.nan)
retvars = retvars.astype('float')

retvolvars = retvars.copy()
retvolvars.index = pd.to_datetime(retvolvars.index).to_period('M')

retvol = maxretvars.groupby(by=maxretvars.index).std()
retvol.index = retvol.index.asfreq('D')

retvol.astype(float).to_csv(pathto_1_200 + '73.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 6 baspread ( Bid-Ask Spread )
'''mean((askhi-bidlo)/((askhi+bidlo)/2)) as baspread
'''
askhivars = pd.read_csv(pathtorawdata + 'dailyaskhi.csv')
askhivars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
askhivars = askhivars.pivot(index='date', columns='PERMNO', values='ASKHI')
askhivars.index = pd.to_datetime(askhivars.index, format='%d/%m/%Y')
askhivars = askhivars.sort_index()

bidlovars = pd.read_csv(pathtorawdata + 'dailybidlo.csv')
bidlovars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
bidlovars = bidlovars.pivot(index='date', columns='PERMNO', values='BIDLO')
bidlovars.index = pd.to_datetime(bidlovars.index, format='%d/%m/%Y')
bidlovars = bidlovars.sort_index()

baspreadvars = ( askhivars - bidlovars ) / ( ( askhivars + bidlovars ) / 2 ).replace(0,np.nan)
baspreadvars.index = pd.to_datetime(baspreadvars.index).to_period('M')

baspread = baspreadvars.groupby(by=baspreadvars.index).mean()

baspread.index = baspread.index.asfreq('D')

baspread.astype(float).to_csv(pathto_1_200 + '6.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 87 std_dolvol ( Volatility of Liquidity ( Dollar Trading Volume ) ) #81
'''std(log(abs(prc*vol))) as std_dolvol
'''
prcvars = pd.read_csv(pathtorawdata + 'dailyprc.csv')
prcvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
prcvars = prcvars.pivot(index='date', columns='PERMNO', values='PRC')
prcvars.index = pd.to_datetime(prcvars.index, format='%d/%m/%Y')
prcvars = prcvars.sort_index()

volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

dolvolvars = np.log(( prcvars * volvars ).abs())
dolvolvars.index = pd.to_datetime(dolvolvars.index).to_period('M')

std_dolvol = dolvolvars.groupby(by=dolvolvars.index).std()

std_dolvol.index = std_dolvol.index.asfreq('D')

std_dolvol.astype(float).to_csv(pathto_1_200 + '87.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 88 std_turn ( Volatility of Liquidity ( Share Turnover ) ) #82
'''std(vol/shrout) as std_turn
'''
volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

shroutvars = pd.read_csv(pathtorawdata + 'dailyshrout.csv')
shroutvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
shroutvars = shroutvars.pivot(index='date', columns='PERMNO', values='SHROUT')
shroutvars.index = pd.to_datetime(shroutvars.index, format='%d/%m/%Y')
shroutvars = shroutvars.sort_index()

std_turnvars = volvars / shroutvars.replace(0,np.nan)
std_turnvars.index = pd.to_datetime(std_turnvars.index).to_period('M')

std_turn = std_turnvars.groupby(by=std_turnvars.index).std()

std_turn.index = std_turn.index.asfreq('D')

std_turn.astype(float).to_csv(pathto_1_200 + '88.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 40 ill ( İlliquidity ) #83
'''mean(abs(ret)/(abs(prc)*vol)) as ill
'''
retvars = pd.read_csv(pathtorawdata + 'dailyret.csv')
retvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
retvars = retvars.pivot(index='date', columns='PERMNO', values='RET')
retvars.index = pd.to_datetime(retvars.index, format='%d/%m/%Y')
retvars = retvars.sort_index()

retvars = retvars.replace('C', np.nan)
retvars = retvars.replace('B', np.nan)
retvars = retvars.astype('float')

prcvars = pd.read_csv(pathtorawdata + 'dailyprc.csv')
prcvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
prcvars = prcvars.pivot(index='date', columns='PERMNO', values='PRC')
prcvars.index = pd.to_datetime(prcvars.index, format='%d/%m/%Y')
prcvars = prcvars.sort_index()

volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

ill = retvars.abs() / ( prcvars.abs() * volvars ).replace(0,np.nan)
ill.index = pd.to_datetime(ill.index).to_period('M')

ill = ill.groupby(by=ill.index).mean()

ill.index = ill.index.asfreq('D')

ill.astype(float).to_csv(pathto_1_200 + '40.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 34
'''sum(vol=0) as countzero
'''
volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

countzerovars = volvars.copy()

countzerovars[countzerovars!=0] = np.nan
countzerovars[countzerovars==0] = 1

countzerovars.index = pd.to_datetime(countzerovars.index).to_period('M')
countzero = countzerovars.groupby(by=countzerovars.index).sum()

countzero.index = countzero.index.asfreq('D')

countzero.astype(float).to_csv(pathtoaddvars_1_200 + '34.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # Additional Variable 35
'''n(permno) as ndays
'''
ndaysvars = pd.read_csv((pathtorawdata + 'daily_crsp.csv'), usecols=['PERMNO','date'])
ndaysvars = dum.copy()
ndaysvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
ndaysvars.set_index('date', inplace=True)
ndaysvars.index = pd.to_datetime(ndaysvars.index, format='%d/%m/%Y').to_period('M')

ndaysvars2 = ndaysvars.groupby(by=[ndaysvars.index,ndaysvars.iloc[:,0]]).count()
ndaysvars2 =ndaysvars2.unstack()

ndays = ndaysvars2.copy()
ndays.index = ndays.index.asfreq('D')

ndays.columns = list(ndaysvars['PERMNO'].unique())

ndays.astype(float).to_csv(pathtoaddvars_1_200 + '35.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 90 turn ( Share Turnover ) #84
'''sum(vol/shrout) as turn
'''
volvars = pd.read_csv(pathtorawdata + 'dailyvol.csv')
volvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
volvars = volvars.pivot(index='date', columns='PERMNO', values='VOL')
volvars.index = pd.to_datetime(volvars.index, format='%d/%m/%Y')
volvars = volvars.sort_index()

shroutvars = pd.read_csv(pathtorawdata + 'dailyshrout.csv')
shroutvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
shroutvars = shroutvars.pivot(index='date', columns='PERMNO', values='SHROUT')
shroutvars.index = pd.to_datetime(shroutvars.index, format='%d/%m/%Y')
shroutvars = shroutvars.sort_index()

turnvars = ( volvars / shroutvars.replace(0,np.nan) )
turnvars.index = pd.to_datetime(turnvars.index).to_period('M')

turn = turnvars.groupby(by=turnvars.index).sum().replace(0,np.nan)

turn.index = turn.index.asfreq('D')

turn.astype(float).to_csv(pathto_1_200 + '90.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 94 zerotrade ( Zero Trading Days ) #85
'''zerotrade=(countzero+((1/turn)/480000))*21/ndays
Use Additional Variable 34 countzero 2166
Use Additional Variable 35 ndays 2183
Use Variable #90 2200
'''
countzero = pd.read_csv(pathtoaddvars_1_200 + '34.csv')
countzero.set_index('date', inplace=True)
countzero.index = pd.to_datetime(countzero.index).to_period('D')
countzero.columns = countzero.columns.map(int)

ndays = pd.read_csv(pathtoaddvars_1_200 + '35.csv')
ndays.set_index('date', inplace=True)
ndays.index = pd.to_datetime(ndays.index).to_period('D')
ndays.columns = ndays.columns.map(int)

turn = pd.read_csv(pathto_1_200 + '90.csv')
turn.set_index('date', inplace=True)
turn.index = pd.to_datetime(turn.index).to_period('D')
turn.columns = turn.columns.map(int)

zerotrade = ( countzero + ( ( 1 / turn.replace(0,np.nan) ) / 480000) ) * 21 / ndays.replace(0,np.nan)

zerotrade.astype(float).to_csv(pathto_1_200 + '94.csv')

#----------------------------------------------------------------------------------------------------------------------------------
    # 7 beta ( Beta ) #86
    # 39 idiovol ( Idiosyncratic Return Volatility ) # 87
'''
'''
'''Raw Data Preperation Process for Beta Estimation'''
'''
# =============================================================================
retvars = pd.read_csv(pathtorawdata + 'dailyret.csv')
retvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
retvars = retvars.pivot(index='date', columns='PERMNO', values='RET')
retvars.index = pd.to_datetime(retvars.index, format='%d/%m/%Y')
retvars = retvars.sort_index()

retvars = retvars.replace('C', np.nan)
retvars = retvars.replace('B', np.nan)
retvars = retvars.astype('float')

prcvars = pd.read_csv(pathtorawdata + 'dailyprc.csv')
prcvars.drop_duplicates(subset=['PERMNO', 'date'], inplace=True)
prcvars = prcvars.pivot(index='date', columns='PERMNO', values='PRC')
prcvars.index = pd.to_datetime(prcvars.index, format='%d/%m/%Y')
prcvars = prcvars.sort_index()

retvars[prcvars<=0] = np.nan

mktrfvars = pd.read_csv(pathtorawdata + 'dailyff.csv')
mktrfvars.set_index('date', inplace=True)
mktrfvars.index = pd.to_datetime(mktrfvars.index, format='%d/%m/%Y')

returnsvars = pd.concat([retvars, mktrfvars], axis=1)
returnsvars.astype(float).to_csv(pathtorawdata + 'variablesforbetaestimation.csv')
# =============================================================================
'''

returnsvars = pd.read_csv(pathtorawdata + 'variablesforbetaestimation.csv')
returnsvars.set_index('date', inplace=True)
ys = returnsvars.iloc[:,:-5]
x = returnsvars['mktrf']
rf = returnsvars['rf']

ys = pd.DataFrame(np.subtract(ys.values,rf.values.reshape(-1,1)), index = ys.index, columns=ys.columns)

''' Variable 1 for beta prices - rf'''
wkret = np.exp((np.log(ys+1)).rolling(5, min_periods=4).sum()) - 1
wkret.index = pd.to_datetime(wkret.index).to_period('W')
wkret = wkret.groupby(by=wkret.index).last()
wkret.index = wkret.index.asfreq('D')
wkret = wkret.loc[:'2018-12-31', :]

#nadum = wkret.isna()

ret = wkret.rolling( (156), min_periods=52 ).mean()

''' Variable 2 for beta market - rf '''
xwkret = np.exp((np.log(x+1)).rolling(5).sum()) -1      
xwkret.index = pd.to_datetime(xwkret.index).to_period('W')
xwkret = xwkret.groupby(by=xwkret.index).last()
xwkret.index = xwkret.index.asfreq('D')
xwkret = pd.DataFrame(xwkret).loc[:'2018-12-31',:]

mktrf = xwkret.rolling(156).mean()
        
''' variable 3 for beta Market * Returns'''
mktrf_retvars = pd.DataFrame( ( wkret.values * xwkret.values.reshape(-1,1) ), index=wkret.index, columns=wkret.columns)
mktrf_ret = mktrf_retvars.rolling( (156), min_periods=52 ).mean()

''' variable 4 for beta Market * Market'''
mktrf2vars = xwkret**2
mktrf2 = mktrf2vars.rolling( (156), min_periods=52 ).mean()

''' variable 5 for beta Returns * Returns'''
ret2vars = wkret**2
ret2 = ret2vars.rolling( (156), min_periods=52 ).mean()

''' variable 6 for beta nandummy '''
nadumvars = pd.DataFrame(np.ones(shape=wkret.shape), index=wkret.index, columns=wkret.columns)
nadumvars[wkret.isna()] = 0
nadum = nadumvars.rolling(156).apply(np.sum)
nadum[nadum<52] = np.nan
nadum[nadum.notnull()] = 1

''' nandummy for missing return days'''
nadumvars2 = pd.DataFrame(np.ones(shape=wkret.shape), index=wkret.index, columns=wkret.columns)
nadumvars2[wkret.isna()] = np.nan
nadum2 = nadumvars2.copy()

''' Beta '''
beta = pd.DataFrame( ( ( mktrf_ret.values - ( ret.values * mktrf.values.reshape(-1,1) ) ) / ( mktrf2.values.reshape(-1,1) - ( mktrf.values.reshape(-1,1) ** 2 ) ) ), index = ret.index, columns=ret.columns )
#beta = beta * nadum
beta = beta * nadum
beta = beta * nadum2

alpha = ret - ( beta * mktrf_ret )

idvol = ret2 - ( alpha**2 ) - ( 2 * alpha * beta * mktrf_ret )

betasq = beta ** 2 


beta.index = (beta.index.to_timestamp()).to_period('M')
beta = beta.groupby(by=beta.index).last()
beta.index = beta.index.asfreq('D')

idvol.index = (idvol.index.to_timestamp()).to_period('M')
idvol = idvol.groupby(by=idvol.index).last()
idvol.index = idvol.index.asfreq('D')

betasq.index = (betasq.index.to_timestamp()).to_period('M')
betasq = betasq.groupby(by=betasq.index).last()
betasq.index = betasq.index.asfreq('D')

beta.astype(float).to_csv(pathto_1_200 + '7.csv')
idvol.astype(float).to_csv(pathto_1_200 + '39.csv')
betasq.astype(float).to_csv(pathto_1_200 + '8.csv')


