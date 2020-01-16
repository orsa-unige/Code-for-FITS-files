#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[1]:


#!/usr/bin/env python3


# In[2]:


# -*- coding: utf-8 -*-


# In[3]:


from astropy.io import fits 
import os
import glob
import click
import path


# In[4]:


cammino = click.prompt('Please provide path: ', type=click.Path(exists=True, resolve_path=True))
print(cammino)
path = os.chdir(cammino)


# In[5]:


x = glob.glob("*.fits")
for files in x:
    hdulist = fits.open(files) 
    hdulist.info()
    hdu = hdulist[0]
    hdu.data.shape


# In[13]:


y = input('Please provide header: ')
hdu.header[y]
quote = y
z = input('Please provide header you want to subsitute with: ')
hdu.header[z]
w = hdu.header[z]
for files in x:
    if(quote.find(y)):
       print("Header already exists")
    else:
        print("Header not found - creating now")
        fits.setval(files, y, value='some value') 
        with fits.open(files, 'update') as f:
            for hdu in f:
                hdu.header[y] = w
                
print("Done!")
