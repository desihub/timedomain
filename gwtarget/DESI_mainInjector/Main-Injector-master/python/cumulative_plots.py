#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys, getopt, os

def usage():
    print 'Usage: cumulative_plots.py [-d <path_to_files>] [-p <path_for_plots> ] -e <event_id> -f <file 1 file 2 ... file n>'

try:
    opts, args = getopt.getopt(sys.argv[1:],"hd:p:e:f:",["files=","event="])
except getopt.GetoptError:
    usage()
    sys.exit(2)

if len(args) > 0:
    usage()
    sys.exit(2)

path = '.'
ppath = '.'
this_event = None
this_file = None
events = [ ]
probs = []

for opt, arg in opts:
    if opt == '-h':
        usage()
        sys.exit()
    elif opt == ('-d'):
        path = arg
        for filename in os.listdir(path):

            events.append(filename.split("-")[0])
            probs.append( np.genfromtxt(path + '/' + filename, usecols=2) )

    elif opt == '-p':
        ppath = arg
    elif opt in ("-e","--event"):
        this_event = arg
    elif opt in ("-f","--files"):
        if len(arg.split()) != 1:
            usage()
            sys.exit(2)
        this_file = arg

events.append(this_file.split("/")[-1].split("-")[0])

try:
    probs.append(np.genfromtxt(this_file, usecols=3))
except IOError:
    print "warning: file for this event not found."

i = dict(zip(events,range(len(events))))

plt.figure(figsize=(8.5*1.618,8.5))
sims_labeled = False

for event in events:
#    print(event, i[event], len(probs))
    cumprobs = np.sort(100*probs[i[event]])
    cumprobs = np.cumsum(cumprobs[::-1])
    if i[event] == i[this_event]:
        plt.plot(cumprobs,color='blue',label='Event #: '+event,linewidth=3.0)
    elif sims_labeled:
        plt.plot(cumprobs,color='black',alpha=0.5,label='_')
    else:        
        plt.plot(cumprobs,color='black',alpha=0.5,label='Simulations')
        sims_labeled = True
    
plt.ylabel('Probability (%)')
plt.xlabel('# of hexes')
plt.title('Cumulative Probability Distribution')
plt.grid(True)
plt.legend()
plt.savefig(ppath+'/'+this_event+'-and-sim-cumprobs.png')
plt.close()



