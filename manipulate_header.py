#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[1]:


#!/usr/bin/env python3


# In[2]:


# -*- coding: utf-8 -*-


# In[41]:


from astropy.io import fits 
import glob
import sys


# In[42]:





# In[56]:


def main():
    for argument in sys.argv:
        print(argument)
        key = sys.argv[1]
        valore = sys.argv[2]
        folder = glob.glob("./**/*.fits", recursive = True)
        for files in folder:
            print(files)
            hdu = fits.PrimaryHDU()
            hdulist = fits.open(files) 
            x = (repr(fits.getheader(files, 0)))
            if key in x:
                print("Header already exists")
            else:
                print("Header not found - creating now")
                fits.setval(files, key, value ='some value')
                with fits.open(files, 'update') as f:
                    for hdu in f:
                        hdu.header[key] = valore

if __name__ == "__main__":
    main()


# In[57]:





# In[ ]:




