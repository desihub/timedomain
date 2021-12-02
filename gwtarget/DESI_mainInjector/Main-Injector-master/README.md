# Main-Injector

As of June 2018, run as user gw in the gw account. 

We'd like to be able to run it in test mode as user XXX (say alenon, annis, brout, marcelle). This entails having the directory structure in the gw account being listed in the maininjector yaml configuration file.

```unix
% git clone git@github.com:SSantosLab/Main-Injector.git
% cd Main-Injector
```

You'll need to run the setup script each time you login
```
% source SOURCEME
```

Then you can run the recycler (this is automatically in test mode)
```
% python recycler.py
```


