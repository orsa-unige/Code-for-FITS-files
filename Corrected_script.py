#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[1]:


#!/usr/bin/env python3


# In[2]:


# -*- coding: utf-8 -*-


# In[7]:


from astropy.io import fits
import os


# In[8]:


path = os.chdir("Dottorato/Test_CCD/03_12_19")


# In[10]:


file1 = "NUNKI.2019-12-03T13_41_21.731Z.fits"
file2 = "NUNKI.2019-12-03T13_46_46.171Z.fits"
file3 = "NUNKI.2019-12-03T13_47_55.836Z.fits"
file4 = "NUNKI.2019-12-03T13_48_08.293Z.fits"
file5 = "NUNKI.2019-12-03T13_48_21.361Z.fits"
file6 = "NUNKI.2019-12-03T13_48_33.651Z.fits"
file7 = "NUNKI.2019-12-03T13_48_47.136Z.fits"
file8 = "NUNKI.2019-12-03T13_48_59.083Z.fits"
file9 = "NUNKI.2019-12-03T13_49_26.879Z.fits"
file10 = "NUNKI.2019-12-03T13_49_39.817Z.fits"
file11 = "NUNKI.2019-12-03T13_50_09.802Z.fits"
file12 = "NUNKI.2019-12-03T13_50_19.248Z.fits"
file13 = "NUNKI.2019-12-03T13_50_34.918Z.fits"
file14 = "NUNKI.2019-12-03T13_50_45.214Z.fits"
file15 = "NUNKI.2019-12-03T13_50_58.754Z.fits"
file16 = "NUNKI.2019-12-03T13_51_12.984Z.fits"
file17 = "NUNKI.2019-12-03T13_51_22.551Z.fits"
file18 = "NUNKI.2019-12-03T13_51_35.746Z.fits"
file19 = "NUNKI.2019-12-03T13_51_48.882Z.fits"
file20 = "NUNKI.2019-12-03T13_52_01.703Z.fits"
file21 = "NUNKI.2019-12-03T13_52_14.214Z.fits"
file22 = "NUNKI.2019-12-03T13_52_27.409Z.fits"
file23 = "NUNKI.2019-12-03T13_52_40.347Z.fits"
file24 = "NUNKI.2019-12-03T13_52_54.941Z.fits"
file25 = "NUNKI.2019-12-03T13_53_09.045Z.fits"
file26 = "NUNKI.2019-12-03T13_53_22.550Z.fits"
file27 = "NUNKI.2019-12-03T13_53_35.242Z.fits"
file28 = "NUNKI.2019-12-03T13_53_48.035Z.fits"
file29 = "NUNKI.2019-12-03T13_54_00.975Z.fits"
file30 = "NUNKI.2019-12-03T13_54_15.075Z.fits"
file31 = "NUNKI.2019-12-03T13_54_26.822Z.fits"
file32 = "NUNKI.2019-12-03T13_54_40.723Z.fits"
file33 = "NUNKI.2019-12-03T13_54_53.801Z.fits"
file34 = "NUNKI.2019-12-03T13_55_06.322Z.fits"
file35 = "NUNKI.2019-12-03T13_55_21.105Z.fits"
file36 = "NUNKI.2019-12-03T13_55_32.815Z.fits"
file37 = "NUNKI.2019-12-03T14_01_08.393Z.fits"
file38 = "NUNKI.2019-12-03T14_01_16.140Z.fits"
file39 = "NUNKI.2019-12-03T14_01_23.875Z.fits"
file40 = "NUNKI.2019-12-03T14_01_31.676Z.fits"
file41 = "NUNKI.2019-12-03T14_01_39.463Z.fits"
file42 = "NUNKI.2019-12-03T14_01_47.215Z.fits"
file43 = "NUNKI.2019-12-03T14_01_55.031Z.fits"
file44 = "NUNKI.2019-12-03T14_02_02.769Z.fits"
file45 = "NUNKI.2019-12-03T14_02_10.538Z.fits"
file46 = "NUNKI.2019-12-03T14_02_18.344Z.fits"
file47 = "NUNKI.2019-12-03T14_02_26.116Z.fits"
file48 = "NUNKI.2019-12-03T14_02_49.381Z.fits"
file49 = "NUNKI.2019-12-03T14_02_57.136Z.fits"
file50 = "NUNKI.2019-12-03T14_03_05.023Z.fits"
file51 = "NUNKI.2019-12-03T14_03_12.764Z.fits"
file52 = "NUNKI.2019-12-03T14_03_20.561Z.fits"
file53 = "NUNKI.2019-12-03T14_03_28.368Z.fits"
file54 = "NUNKI.2019-12-03T14_03_36.097Z.fits"
file55 = "NUNKI.2019-12-03T14_03_43.903Z.fits"
file56 = "NUNKI.2019-12-03T14_03_51.676Z.fits"
file57 = "NUNKI.2019-12-03T14_03_59.416Z.fits"
file58 = "NUNKI.2019-12-03T14_04_07.213Z.fits"
file59 = "NUNKI.2019-12-03T14_04_14.926Z.fits"
file60 = "NUNKI.2019-12-03T14_04_22.714Z.fits"
file61 = "NUNKI.2019-12-03T14_04_30.500Z.fits"
file62 = "NUNKI.2019-12-03T14_04_38.226Z.fits"
file63 = "NUNKI.2019-12-03T14_04_46.016Z.fits"
file64 = "NUNKI.2019-12-03T14_04_53.744Z.fits"
file65 = "NUNKI.2019-12-03T14_12_38.760Z.fits"
file66 = "NUNKI.2019-12-03T14_19_43.488Z.fits"
file67 = "NUNKI.2019-12-03T14_22_31.044Z.fits"
file68 = "NUNKI.2019-12-03T14_23_28.241Z.fits"
file69 = "NUNKI.2019-12-03T14_25_53.245Z.fits"
file70 = "NUNKI.2019-12-03T14_28_26.703Z.fits"
file71 = "NUNKI.2019-12-03T15_04_15.567Z.fits"
file72 = "NUNKI.2019-12-03T15_09_47.107Z.fits"
file73 = "NUNKI.2019-12-03T15_15_28.207Z.fits"
file74 = "NUNKI.2019-12-03T15_23_52.141Z.fits"
file75 = "NUNKI.2019-12-03T15_34_18.057Z.fits"
file76 = "NUNKI.2019-12-03T15_39_52.687Z.fits"


# In[11]:


hdulist = fits.open(file1) 


# In[12]:


hdulist.info()
filename = file1


# In[15]:


hdu = hdulist[0]


# In[16]:


hdu.data.shape


# In[17]:


hdu.header


# In[18]:


fits.setval(file1, 'IMAGETYP', value='some value')


# In[19]:


with fits.open(file1, 'update') as f:
    for hdu in f:
        hdu.header['IMAGETYP'] = 'Dark'

       
print("After modifications:")
print()
print("Extension 0:")
print(repr(fits.getheader(file1, 0)))
print()


