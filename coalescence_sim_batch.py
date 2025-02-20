#!/usr/bin/env python
# coding: utf-8

# In[ ]:


FST = 0.1


# In[ ]:


for event in ("split", "migration"):
    if event == "split":
        for subtype in ("star",):
            for d in (2,4,8,16):
                print("Running %s model, %s with %s demes" % (event, subtype, d))
                get_ipython().run_line_magic('run', 'coalescence_sim.ipynb')
        for subtype in ("balanced", "caterpillar"):
            for d in (4,8,16):
                print("Running %s model, %s with %s demes" % (event, subtype, d))
                get_ipython().run_line_magic('run', 'coalescence_sim.ipynb')
    if event == "migration":
        for subtype in ("island",):
            for d in (2,4,8,16):
                print("Running %s model, %s with %s demes" % (event, subtype, d))
                get_ipython().run_line_magic('run', 'coalescence_sim.ipynb')
        for subtype in ("circular",):
            for d in (4,8,16):
                print("Running %s model, %s with %s demes" % (event, subtype, d))
                get_ipython().run_line_magic('run', 'coalescence_sim.ipynb')


# In[ ]:




