import numpy as np
from email.mime.text import MIMEText
import smtplib
import requests

def postImagetoSLACK(impath,official=False):
    if not official:
        return
    post = {}
    post["value1"] = open(impath,'rb').read()
    post["value2"] = '' 
    r = requests.post('https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh',data=post)

def postToSLACK(subject,text,official=False,atchannel=False):
    print(official)
    if not official:
        return
    if atchannel:
        subject = '@channel ' + subject
        

    post = {}
    post["value1"] = subject
    post["value2"] = text
    requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=post)

def send(subject,text,official=False,atchannel=False):
    #import MIalerts as config
    #del config
    #import MIalerts as config
    #asdf
    if not official:
        return

    msg = MIMEText(text)
    
    try:
        people = np.genfromtxt('/home/s1/desgw/PHONE_AND_EMAIL_LIST.TXT', dtype=[('name','S50'),('email','S50'),('phone','S50')], delimiter=",",skip_header=1)
    except:
        report = {}
        report["value1"] = '@channel WARNING, THE PHONE AND EMAIL LIST CRASHED... probably a formatting issue \n'
        report["value2"] = 'please fix asap.'
        requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=report)

    print people['name']

    me = 'automated-desGW@fnal.gov'
    you = []
    if official:
        for email in people['email']:
            #if not p[0] == '#':
            you.append(email.replace('\t','').replace(' ',''))
    #else:
    #    you = ['djbrout@gmail.com']

    if official:
        for phone in people['phone']:
            if not phone == 'none':
                #if not p[0] == '#':
                you.append(phone.replace('\t','').replace(' ',''))

        report = {}
        if atchannel:
            subject = '@channel '+subject
        report["value1"] = subject
        report["value2"] = text
        requests.post("https://maker.ifttt.com/trigger/desgwtrigger/with/key/dmCI74bFTHzz2uCOIevguh", data=report)            

    #else:
    #    you.append('2153008763@mms.att.net')

    for y in you:
        try:
            msg['Subject'] = subject
            msg['From'] = me

            #msg['From'] = me
            msg['To'] = y

            s = smtplib.SMTP('localhost')
            s.sendmail(me, y, msg.as_string())
            s.quit()
        except:
            print 'could not send alert! to ', y

        

#send('This is just a test #6 DO NOT REPLY VIA EMAIL','LOOK AT SLACK PAGE, cross your fingers, should be instant...')
