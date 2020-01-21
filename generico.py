#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[1]:


#!/usr/bin/env python3


# In[2]:


# -*- coding: utf-8 -*-


# In[3]:


from astropy.io import fits 
import os
import click
import path


# In[4]:


cammino = click.prompt('Please provide path: ', type=click.Path(exists=True, resolve_path=True))
print(cammino)
path = os.chdir(cammino)
file = input('Please provide file.fits: ')
print(file)


# In[5]:


hdu = fits.PrimaryHDU()
hdulist = fits.open(file) 
hdulist.info()


# In[15]:


y = input("Please provide header: ")
mylist = list(hdu.header)
print(y)
if y in mylist:
    print("Header already exists")
else:
    print("Header not found - creating now")
    fits.setval(file, y, value='some value')
    z = input("Please provide header's value: ")
    with fits.open(file, 'update') as f:
        for hdu in f:
            hdu.header[y] = z
print(hdu.header)

