import os
import numpy as np

import datetime
from datetime import datetime as dt
import pytz

def mjd_to_datetime(mjd):
    mjd_epoch = datetime.datetime(1858, 11, 17, tzinfo=pytz.utc)
    d = mjd_epoch + datetime.timedelta(mjd)
    return d

def nowUTC():
    now = dt.utcnow()
    return now

def makeNewPage(master_dir,outfilename,trigger_id,event_paramfile,mapfolder,processing_param_file=None,real_or_sim='real',secondtimearound=False):
    #print 'inside triggerpages'
    #print event_paramfile
    #print 'making new page'*100
    #print outfilename
    #print trigger_id
    #print mapfolder
    #print secondtimearound

    try:
        event_params = np.load(event_paramfile)
    except:
        print('could not load event_paramfile',event_paramfile)
        return
    d = mjd_to_datetime(float(str(event_params['MJD'])))
    n = nowUTC()
    try:
        processing = np.load(processing_param_file)
    except:
        processing = None

    #if secondtimearound:
    #tmapfolder = '.'
    #else:
    tmapfolder = mapfolder

    candidates = {'id':[],'ra':[],'dec':[],'peakmag':[],'peakmjd':[],'mlscore':[]}
    # try:
    #     candidate_param_files = os.listdir('DES_GW_Website/Triggers/'+trigger_id+'/candidate_param_files/')
    #
    #     if len(candidate_param_files) > 0:
    #         for cpf in candidate_param_files:
    #             cp = np.load(cpf)
    #             candidates['id'].append(cp['id'])
    #             candidates['ra'].append(cp['ra'])
    #             candidates['dec'].append(cp['dec'])
    #             candidates['peakmag'].append(cp['peakmag'])
    #             candidates['peakmjd'].append(cp['peakmjd'])
    #             candidates['mlscore'].append(cp['peakmlscore'])
    # except OSError:
    #     candidate_param_files = []
    html = \
        '<!DOCTYPE html> \
        <html lang="en" style="height:100%;">\
        <head> \
        <meta charset="utf-8"> \
        <title>DES GW</title>\
        <meta name="viewport" content="width=device-width, initial-scale=1.0"> \
        <meta name="keywords" content="pinegrow, blocks, bootstrap" />\
        <meta name="description" content="DESGW" />\
        <link rel="shortcut icon" href="ico/favicon.png"> \
        <!-- Core CSS -->\
        <link href="/desgw/assets/bootstrap/css/bootstrap.min.css" rel="stylesheet"> \
        <link href="/desgw/assets/css/font-awesome.min.css" rel="stylesheet">\
        <link href="http://fonts.googleapis.com/css?family=Open+Sans:300italic,400italic,600italic,700italic,400,300,600,700" rel="stylesheet">\
        <link href="http://fonts.googleapis.com/css?family=Lato:300,400,700,300italic,400italic,700italic" rel="stylesheet">\
        <!-- Style Library -->         \
        <link href="/desgw/assets/css/style-library-1.css" rel="stylesheet">\
        <link href="/desgw/assets/css/plugins.css" rel="stylesheet">\
        <link href="/desgw/assets/css/blocks.css" rel="stylesheet">\
        <link href="/desgw/assets/css/custom.css" rel="stylesheet">\
        <!-- HTML5 shim, for IE6-8 support of HTML5 elements. All other JS at the end of file. -->         \
        <!--[if lt IE 9]>\
        <script src="/desgw/assets/js/html5shiv.js"></script>\
        <script src="/desgw/assets/js/respond.min.js"></script>\
            <script type="text/javascript" src="sortable.js"></script>\
            <![endif]-->\
                <script type="text/javascript">\
            function ChangeColor(tableRow, highLight)\
            {\
            if (highLight)\
            {\
              tableRow.style.backgroundColor = "#7FFFD4";\
            }\
            else\
            {\
              tableRow.style.backgroundColor = "white";\
            }\
          }\
          function DoNav(theUrl)\
          {\
          document.location.href = theUrl;\
          }\
          </script>\
            </head>' + \
            '<body data-spy="scroll" data-target="nav"> \
                <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>\
                <script type="text/javascript" src="/desgw/assets/js/smoothscroll.js"></script>\
                <script type="text/javascript" src="/desgw/assets/js/breedjs.min.js"></script>\
                <header id="header-2" class="soft-scroll header-2">\
                    <nav class="main-nav navbar navbar-default navbar-fixed-top">\
                        <div class="container">\
                            <!-- Brand and toggle get grouped for better mobile display -->\
                            <div class="navbar-header">\
                                <a href="/desgw">\
                                    <button type="button" style="outline: none;" class="btn btn-default">Home</button>\
                                </a>\
                                <button type="button" style="outline: none;" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar-collapse">\
                                    <span class="sr-only">Toggle navigation</span>\
                                    <span class="icon-bar"></span>\
                                    <span class="icon-bar"></span>\
                                    <span class="icon-bar"></span>\
                                </button>\
                                <a href="#"></a>\
                            </div>\
                            <!-- Collect the nav links, forms, and other content for toggling -->\
                            <div class="collapse navbar-collapse" id="navbar-collapse">\
                                <ul class="nav navbar-nav navbar-right">\
                                    <li class="nav-item active">\
                                        <a href="#Trigger" class="smoothScroll">'+str(trigger_id)+'</a>\
                                    </li>\
                                    <li class="nav-item">\
                                        <a href="#LIGO" class="smoothScroll">LIGO Data</a>\
                                    </li>\
                                    <li class="nav-item">\
                                        <a href="#DES" class="smoothScroll">DES Data</a>\
                                    </li>\
                                    <li class="nav-item">\
                                        <a href="#Plots" class="smoothScroll">Plots</a>\
                                    </li>\
                                    <li class="nav-item">\
                                        <a href="../../post-processing-all/Master-master.html" class="smoothScroll">Processing</a>\
                                    </li>\
                                    <li class="nav-item">\
                                        <a href="../../imalive.html" class="smoothScroll">Alive?</a>\
                                    </li>\
                                </ul>\
                            </div>\
                            <!-- /.navbar-collapse -->\
                        </div>\
                        <!-- /.container-fluid -->\
                    </nav>\
                </header>\
                <section id="promo-3" class="content-block promo-3 min-height-600px bg-deepocean">\
                    <div class="container text-center">\
                        <div class="gallery-item" style="width: 50%; left: 50%; margin-right: -50%; transform: translate(50%, 0%)">\
                         <div class="gallery-thumb" >\
                            <img name="Trigger" src="'+tmapfolder+'/'+str(trigger_id)+'_animate.gif" onerror="this.src=\'/desgw/image_placeholder.jpg\'" \
                            class="img-responsive" alt="2nd gallery Thumb" width="50%">\
                            <div class="image-overlay"></div>\
                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_animate.gif" onerror="this.href=\'/desgw/image_placeholder.jpg\'" class="gallery-zoom"><i class="fa fa-eye"></i></a>\
                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_animate.gif" onerror="this.href=\'/desgw/image_placeholder.jpg\'" download class="gallery-link"><i class="fa fa-link"></i></a>\
                          </div>\
                        </div>\
                        <h1>Trigger '+str(trigger_id)+'</h1>\
                        <h4 style="color: white">'+str(d.strftime('%H:%M:%S \t %b %d, %Y UTC'))+'</h4> <h4 style="color: white">FAR: '+str(event_params['FAR'])+'</h4>\
                        <h2></h2><h2>'

    maps = []
    mapsubdirs = []
    print('TMAPFOLDER',tmapfolder,'/',trigger_id,'-observingPlot.gif')
    if real_or_sim == 'real':
        subdir = 'real-triggers'
    if real_or_sim == 'sim':
        subdir = 'test-triggers'
    print(master_dir)
    for fff in os.listdir(master_dir):
        if '.processing' in fff:
            maps.append(fff.split('.')[0])
            
    #os.system('scp ./'+subdir+'/'+trigger_id+'/'+map+'/'+trigger_id+'_trigger.html codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + self.trigger_id + '/')
            
    thismap = ''
    for map in maps:
        print(map)
        #raw_input()
        print(list(event_params.keys()))
        if map in str(event_params['filename']):
                html+='<div class="btn-group">\
                             <a ><button\
                            type = "button"\
                            style="width:260px; background-color:white; color:black; outline:none;"\
                            class="btn btn-outline btn-outline-xl outline-light"> '+map+' </button></a>'
                thismap = map
                epf = './'+subdir+'/'+trigger_id+'/'+map+'/'+trigger_id+'_params.npz'
                ofn = './'+subdir+'/'+trigger_id+'/'+trigger_id+'_'+map+'_trigger.html'
                if not secondtimearound:
                    makeNewPage(master_dir,ofn,trigger_id,epf,map,real_or_sim=real_or_sim,secondtimearound=True)
        else:
                html += '<div class="btn-group">\
                    <a href="'+trigger_id+'_'+map+'_trigger.html"><button\
                    type="button"\
                    style="width:260px; background-color:transparent; color:white; outline:none;"\
                    class="btn btn-outline btn-outline-xl outline-light"> '+map+' </button></a>'
                if not secondtimearound:
                    epf = './'+subdir+'/'+trigger_id+'/'+map+'/'+trigger_id+'_params.npz'
                    ofn = './'+subdir+'/'+trigger_id+'/'+trigger_id+'_'+map+'_trigger.html'
                    makeNewPage(master_dir,ofn,trigger_id,epf,map,real_or_sim=real_or_sim,secondtimearound=True)

    html += '</h2>'
    html +='<h2>LIGO Probability of Detection (in our hexes): ' + str(round(float(str(event_params['LIGO_prob'])), 6)) + '</h2>\
                        <h2>Trigger Type: '+str(event_params['boc'])+'</h2>\
                        <h2>Strategy: ' + str(event_params['gethexobstype']) + '</h2>\
    <a href="'+'./'+str(trigger_id)+"_"+thismap+'_JSON.zip" download class="btn btn-outline btn-outline-xl outline-light">Download .json <span class="fa fa-download"></span></a>\
<a href="./'+thismap+'_recycler.log" download class="btn btn-outline btn-outline-xl outline-light">Download .log <span class="fa fa-download"></span></a>\
                    </div>\
                    <!-- /.container -->\
                </section>\
                <h1><a name="LIGO" style="padding-top: 105px;">LIGO DATA</a></h1>\
                <table class="table">\
                    <thead>\
                        <tr>\
                            <th></th>\
                            <th></th>\
                            <th></th>\
                            <th></th>\
                        </tr>\
                    </thead>\
                    <tbody>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>FAR</td>\
                            <td>False Alarm Rate</td>\
                            <td>'+str(event_params['FAR'])+'</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>ETA</td>\
                            <td>Estimated ratio of reduced mass to total mass</td>\
                            <td>'+str(event_params['ETA'])+'</td>\
                        </tr>\
<tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);"><th></th><td>NSProb</td><td>NS Probability</td><td>'+str(event_params['probhasns'])+'</td></tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Chirp Mass</td>\
                            <td>Estimated CBC chirp mass</td>\
                            <td>'+str(event_params['ChirpMass'])+'</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Dist</td>\
                            <td>Estimated distance for CBC event</td>\
                            <td>'+str(round(event_params['codeDistance']))+' Mpc</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>CFreq</td>\
                            <td>Central frequency for chirp event</td>\
                            <td>'+str(event_params['CentralFreq'])+'</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>MJD</td>\
                            <td>Modified Julian Date</td>\
                            <td>'+str(event_params['MJD'])+'</td>\
                        </tr>\
                    </tbody>\
                </table>\
                <h1><a name="DES" style="padding-top: 105px;">DES DATA</a></h1>\
                <table class="table">\
                    <thead>\
                        <tr>\
                            <th></th>\
                            <th></th>\
                            <th></th>\
                            <th></th>\
                        </tr>\
                    </thead>\
                    <tbody>\
                    <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td></td>\
                            <td>Last Processing</td>\
                            <td>'+str(n.strftime('%H:%M:%S \t %b %d, %Y UTC'))+'</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>nHexes</td>\
                            <td>Number of DES Hexes Used</td>\
                            <td>'+str(event_params['nHexes'])+'</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Dist</td>\
                            <td>Model Distance Used</td>\
                            <td>' + str(round(event_params['codeDistance'])) + ' Mpc</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>EProb</td>\
                            <td>Econ Probability</td>\
                            <td>' + str(event_params['econ_prob']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>EArea</td>\
                            <td>Econ Area</td>\
                            <td>' + str(event_params['econ_area']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>NArea</td>\
                            <td>Need Area</td>\
                            <td>' + str(event_params['need_area']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Qual</td>\
                            <td>Quality</td>\
                            <td>' + str(event_params['quality']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>ExpTime</td>\
                            <td>Exposure time in each filter</td>\
                            <td>' + str(event_params['exposure_times']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Filt</td>\
                            <td>Filters used for each exposure</td>\
                            <td>' + str(event_params['exposure_filter']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Hours</td>\
                            <td>Time budget allotted</td>\
                            <td>' + str(event_params['hours']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>Nvisits</td>\
                            <td>Number of visits to area of coverage</td>\
                            <td>' + str(event_params['nvisits']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>NSlots</td>\
                            <td>Number of slots used</td>\
                            <td>' + str(event_params['n_slots']) + '</td>\
                        </tr>\
                        <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                            <th></th>\
                            <td>BSlot</td>\
                            <td>Slot with largest probability</td>\
                            <td>' + str(event_params['best_slot']) + '</td>\
                        </tr>\
                </table>\
                <h1><a name="Plots" style="padding-top: 105px;">Plots</a></h1>\
                <section class="content-block gallery-1 gallery-1-3">\
                    <div class="container">\
                        <!-- /.gallery-filter -->\
                        <div class="row">\
                            <div class="isotope-gallery-container">\
                                <!--<div class="col-sm-6 col-xs-12 gallery-item-wrapper nature outside">\
                                    <!--<div class="gallery-item">\
                                        <!--<div class="gallery-thumb">\
                                            <!--<img src="'+tmapfolder+'/'+str(trigger_id)+'-and-sim-cumprobs.png">\
                                            <!--<div class="image-overlay"></div>\
                                            <!--<a href="'+tmapfolder+'/'+str(trigger_id)+'-and-sim-cumprobs.png" class="gallery-zoom"><i class="fa fa-eye" alt="This is the title"></i></a>\
                                            <!--<a href="'+tmapfolder+'/'+str(trigger_id)+'-and-sim-cumprobs.png" download class="gallery-link" target="_blank"><i class="fa fa-link"></i></a>\
                                        <!--</div>\
                                        <!--<div class="gallery-details">\
                                        <!--    <h5>Cumulative Probability vs Number of Hexes Observed</h5>\
                                        <!--    <p>grey curves from mock data challenge</p>\
                                        <!--</div>\
                                    <!--</div>\
                                <!--</div>-->\
                                <!-- /.gallery-item-wrapper -->\
                                <div class="col-sm-6 col-xs-12 gallery-item-wrapper nature outside"><div class="gallery-item"><div class="gallery-thumb"><img src="'+tmapfolder+'/'+str(trigger_id)+'-slot-probabilities.png" class="img-responsive" alt="2nd gallery Thumb"><div class="image-overlay"></div><a href="'+tmapfolder+'/'+str(trigger_id)+'-slot-probabilities.png" class="gallery-zoom"><i class="fa fa-eye"></i></a><a href="'+tmapfolder+'/'+str(trigger_id)+'-slot-probabilities.png" download class="gallery-link"><i class="fa fa-link"></i></a></div><div class="gallery-details"><h5>DES Probability vs Slot Number</h5><p><br></p></div></div></div><div class="col-sm-6 col-xs-12 gallery-item-wrapper nature outside"><div class="gallery-item"><div class="gallery-thumb"><img src="'+tmapfolder+'/'+str(trigger_id)+'_centered_animate.gif.png" class="img-responsive" alt="2nd gallery Thumb"><div class="image-overlay"></div><a href="'+tmapfolder+'/'+str(trigger_id)+'_centered_animate.gif.png" class="gallery-zoom"><i class="fa fa-eye"></i></a><a href="'+tmapfolder+'/'+str(trigger_id)+'_centered_animate.gif.png" download class="gallery-link"><i class="fa fa-link"></i></a></div><div class="gallery-details"><h5>Zoomed Centered GIF</h5><p><br></p></div></div></div><div class="col-sm-6 col-xs-12 gallery-item-wrapper nature outside">\
                                    <div class="gallery-item">\
                                        <div class="gallery-thumb">\
                                            <img src="'+tmapfolder+'/'+str(trigger_id)+'_LIGO.png" class="img-responsive" alt="2nd gallery Thumb">\
                                            <div class="image-overlay"></div>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_LIGO.png" class="gallery-zoom"><i class="fa fa-eye"></i></a>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_LIGO.png" class="gallery-link"><i class="fa fa-link"></i></a>\
                                        </div>\
                                        <div class="gallery-details">\
                                            <h5>LIGO Probability Contour</h5>\
                                            <p></p>\
                                        </div>\
                                    </div>\
                                </div>\
                                <div class="col-sm-6 col-xs-12 gallery-item-wrapper nature outside">\
                                    <div class="gallery-item">\
                                        <div class="gallery-thumb">\
                                            <img src="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbxLIGO.png" class="img-responsive" alt="2nd gallery Thumb">\
                                            <div class="image-overlay"></div>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbxLIGO.png" class="gallery-zoom"><i class="fa fa-eye"></i></a>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbxLIGO.png" download class="gallery-link"><i class="fa fa-link"></i></a>\
                                        </div>\
                                        <div class="gallery-details">\
                                            <h5>DES Probability X LIGO Probability</h5>\
                                            <p><br></p>\
                                        </div>\
                                    </div>\
                                </div>\
                                <!-- /.gallery-item-wrapper -->\
                                <div class="col-sm-6 col-xs-12 gallery-item-wrapper photography artwork">\
                                    <div class="gallery-item">\
                                        <div class="gallery-thumb">\
                                            <img src="'+tmapfolder+'/'+str(trigger_id)+'_limitingMagMap.png" class="img-responsive" alt="3rd gallery Thumb">\
                                            <div class="image-overlay"></div>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_limitingMagMap.png" class="gallery-zoom"><i class="fa fa-eye"></i></a>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_limitingMagMap.png" download class="gallery-link"><i class="fa fa-link"></i></a>\
                                        </div>\
                                        <div class="gallery-details">\
                                            <h5>DES Limiting Magnitude Map</h5>\
                                            <p><br></p>\
                                        </div>\
                                    </div>\
                                </div>\
                                <!-- /.gallery-item-wrapper -->\
                                <div class="col-sm-6 col-xs-12 gallery-item-wrapper creative nature">\
                                    <div class="gallery-item">\
                                        <div class="gallery-thumb">\
                                            <img src="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbMap.png" class="img-responsive" alt="4th gallery Thumb">\
                                            <div class="image-overlay"></div>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbMap.png" class="gallery-zoom"><i class="fa fa-eye"></i></a>\
                                            <a href="'+tmapfolder+'/'+str(trigger_id)+'_sourceProbMap.png" download class="gallery-link"><i class="fa fa-link"></i></a>\
                                        </div>\
                                        <div class="gallery-details">\
                                            <h5>DES Limiting Mag Map for Source</h5>\
                                            <p>Assumed Source '+str(event_params['codeDistance'])+' Mpc</p>\
                                        </div>\
                                    </div>\
                                </div>\
                                <!-- /.gallery-item-wrapper -->\
                            </div>\
                            <!-- /.isotope-gallery-container -->\
                        </div>\
                        <!-- /.row -->\
                    </div>'





    # processing_table = '\
    #                 <h1><a name="Processing" style="padding-top: 105px;">Processing</a></h1>\
    #                 <table class="table">\
    #                     <thead>\
    #                         <tr>\
    #                             <th></th>\
    #                             <th>GWId</th>\
    #                             <th>Hexnum</th>\
    #                             <th>RA</th>\
    #                             <th>DEC</th>\
    #                             <th>Filter</th>\
    #                             <th>EXPTIME</th>\
    #                             <th>EXPNUM</th>\
    #                             <th>Status</th>\
    #                         </tr>\
    #                     </thead>\
    #                     <tbody>'
    # if not processing is None:
    #
    #     for i,h,mjd,ra,dec,f,et,enum,stat in zip(processing['gwid'],processing['hexnum'],processing['mjd'],
    #                                            processing['ra'],processing['dec'],processing['filt'],
    #                                            processing['exptime'],processing['expnum'],processing['status']):
    #
    #         processing_table += '\
    #                 <tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
    #                     <th></th>\
    #                     <th>'+str(i)+'</th>\
    #                     <th>'+str(h)+'</th>\
    #                     <th>'+str(ra)+'</th>\
    #                     <th>'+str(dec)+'</th>\
    #                     <th>'+str(f)+'</th>\
    #                     <th>'+str(et)+'</th>\
    #                     <th>'+str(enum)+'</th>\
    #                     <th>'+str(stat)+'</th>\
    #                 </tr>'
    #
    # processing_table += '</tbody></table>'
    # html += processing_table
    if not processing is None:

        processing_table = '     <script type="text/javascript">\
                                $(function(){\
                                var data = {\
                                hexes: ['
        for h, mjd, ra, dec, f, et, enum, stat in zip(processing['hexnum'], processing['mjd'],
                                                         processing['ra'], processing['dec'], processing['filt'],
                                                         processing['exptime'], processing['expnum'],
                                                         processing['status']):

            processing_table += '{Hexnum: "'+str(h)+'",MJD: "'+str(mjd)+'", RA: "'+str(ra)+'", DEC: "'+str(dec)+'", Filter: "'+str(f)+'", ' \
                        'EXPTIME: "'+str(et)+'", EXPNUM: "'+str(enum)+'", Status: "'+str(stat)+'" },'

        processing_table += ']\
              };\
              breed.run({\
                scope: "hexes",\
                input: data\
              });\
              $("th").click(function(){\
                breed.sort({\
                  scope: "hexes",\
                  selector: $(this).attr("sort")\
                });\
              });\
            });  </script>\
        <h1><a name="Processing" style="padding-top: 105px;">Processing</a></h1>\
        <table class="table" data-pg-collapsed>\
            <thead>\
                <tr >\
                    <th sort="Hexnum" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Hexnum</th>\
                    <th sort="MJD" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >MJD</th>\
                    <th sort="RA" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >RA</th>\
                    <th sort="DEC" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >DEC</th>\
                    <th sort="Filter" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Filter</th>\
                    <th sort="EXPTIME" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >EXPTIME</th>\
                    <th sort="EXPNUM" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >EXPNUM</th>\
                    <th sort="Status" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Status</th>\
                </tr>\
            </thead>\
            <tbody>\
                <tr b-scope="hexes" b-loop="hex in hexes" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">\
                    <td b-sort="Hexnum">{{hex.Hexnum}}</td>\
                    <td b-sort="MJD">{{hex.MJD}}</td>\
                    <td b-sort="RA">{{hex.RA}}</td>\
                    <td b-sort="DEC">{{hex.DEC}}</td>\
                    <td b-sort="Filter">{{hex.Filter}}</td>\
                    <td b-sort="EXPTIME">{{hex.EXPTIME}}</td>\
                    <td b-sort="EXPNUM">{{hex.EXPNUM}}</td>\
                    <td b-sort="Status">{{hex.Status}}</td>\
                </tr>\
            </tbody>\
        </table>'

        html += processing_table

    if len(candidates['id']) > 0:
        candidates_table = '     <script type="text/javascript">\
                                    $(function(){\
                                    var data = {\
                                    candidates: ['
        for i, ra, dec, pf, pfmjd, ml in zip(candidates['id'], candidates['ra'], candidates['dec'],
                                             candidates['peakmag'], candidates['peakmjd'], candidates['mlscore']):
            candidates_table += '{ID: "' + str(h) + '",RA: "' + str(ra) + '", DEC: "' + str(dec) + '", Peak_Mag: "' +\
                                str(pf) + '", Peak_MJD: "' + str(pfmjd) + '", ML_Score: "' + str(ml) + '" },'

        candidates_table += ']\
                  };\
                  breed.run({\
                    scope: "candidates",\
                    input: data\
                  });\
                  $("th").click(function(){\
                    breed.sort({\
                      scope: "candidates",\
                      selector: $(this).attr("sort")\
                    });\
                  });\
                });  </script>\
            <h1><a name="Candidates" style="padding-top: 105px;">Candidates</a></h1>\
            <table class="table" data-pg-collapsed>\
                <thead>\
                    <tr >\
                        <th sort="ID" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">ID</th>\
                        <th sort="RA" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >RA</th>\
                        <th sort="DEC" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >DEC</th>\
                        <th sort="Peak_Mag" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Peak Mag</th>\
                        <th sort="Peak_MJD" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Peak MJD</th>\
                        <th sort="ML_Score" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >ML Score</th>\
                    </tr>\
                </thead>\
                <tbody>\
                    <tr b-scope="hexes" b-loop="cand in candidates" onmouseover="ChangeColor(this, true);" ' \
                            'onmouseout="ChangeColor(this, false);" onclick="DoNav(http://des-ops.fnal.gov:8080/desgw/Candidates/DES{{cand.ID}}.html);" >\
                        <td b-sort="ID">{{cand.ID}}</td>\
                        <td b-sort="RA">{{cand.RA}}</td>\
                        <td b-sort="DEC">{{cand.DEC}}</td>\
                        <td b-sort="Peak_Mag">{{cand.Peak_Mag}}</td>\
                        <td b-sort="Peak_MJD">{{cand.Peak_MJD}}</td>\
                        <td b-sort="ML_Score">{{cand.ML_Score}}</td>\
                    </tr>\
                </tbody>\
            </table>'

        html += candidates_table


    html += '<!-- /.container -->\
            <div class="copyright-bar bg-white">\
            </div >\
        </section>\
        <script type="text/javascript" src="/desgw/assets/js/jquery-1.11.1.min.js"></script>\
        <script type="text/javascript" src="/desgw/assets/js/bootstrap.min.js"></script>\
        <script type="text/javascript" src="/desgw/assets/js/plugins.js"></script>\
        <script src="https://maps.google.com/maps/api/js?sensor=true"></script>\
        <script type="text/javascript" src="/desgw/assets/js/bskit-scripts.js"></script>\
    </body> \
    </html>'


    a = open(outfilename, 'w')
    a.write(html)
    a.close()
    os.system('scp '+outfilename+' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + trigger_id + '/')
    os.system('scp '+outfilename+' codemanager@desweb.fnal.gov:/des_web/www/html/desgw/Triggers/' + trigger_id + '/'+trigger_id+'_trigger.html')

def make_index_page(webpage_dir, real_or_sim='real'):
    if real_or_sim == 'real':
        fff = 'real-trigger_list.txt'
    if real_or_sim == 'sim':
        fff = 'test-trigger_list.txt'
    tt = open(os.path.join(webpage_dir, fff), 'r')
    lines = tt.readlines()
    tt.close()
    triggers = ''
    isFirst = True
    firstTrigger = ''

    indextable = ''

    index_sorting_script='<script type="text/javascript">\
                            var table = $("table");\
                            $("#trigger, #iprob, #far, #cmass, #mjd, #date")\
                            .wrapInner("<span title="sort this column"/>")\
                            .each(function(){\
                            var th = $(this),\
                            thIndex = th.index(),\
                            inverse = false;\
                            th.click(function(){\
                            table.find("td").filter(function(){\
                            return $(this).index() === thIndex;\
                            }).sortElements(function(a, b){\
                            return $.text([a]) > $.text([b]) ?\
                            inverse ? -1 : 1\
                            : inverse ? 1 : -1;\
                            }, function(){\
                            // parentNode is the element we want to move\
                            return this.parentNode;\
                            });\
                            inverse = !inverse;\
                            });\
                            });  </script>'

    indextable += index_sorting_script
    if real_or_sim == 'real':
        indextable += '<h1><a name="Triggers" style="padding-top: 105px;">Real Triggers</a></h1>'
    else:
        indextable += '<h1><a name="Triggers" style="padding-top: 105px;">Mock Triggers</a></h1>'

    #indextable += '<span><h1><a name="Triggers" style="padding-top: 105px;">Im Alive?</a></h1></span>'

    indextable += '<table class="table" data-pg-collapsed>\
        <thead>\
            <tr >\
                <th></th>\
                <th id="trigger" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);">Trigger</th><th id="type" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Type</th>\
                <th id="iprob" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >LIGO Prob</th>\
                <th id="far" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >FAR</th><th id="dist" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Dist.</th>\
                <th id="mjd" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >MJD</th>\
                <th id="date" onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" >Date</th>\
            </tr>\
        </thead>\
        <tbody>'

    for line in reversed(lines):
        trig, outfolder = line.split(' ')
        if isFirst:
            firstTrigger = trig
            isFirst = False
        outfolder = outfolder.strip('\n')
        params = np.load(os.path.join(outfolder, trig + '_params.npz'))

        d = mjd_to_datetime(float(str(params['MJD'])))
        try:
            #params = np.load(os.path.join(outfolder, trig + '_params2.npz'))
            indextable += '<tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" \
                        onclick="DoNav(\'Triggers/'+trig+'/'+trig+'_trigger.html\');">\
                          <td></td>\
                          <td >'+trig+'</td>\
                          <td >'+str(params['gethexobstype'])+'</td><td >'+str(round(float(str(params['LIGO_prob'])), 6))+'</td>\
                          <td >'+str(params['FAR'])+'</td>\
                          <td >'+str(params['codeDistance'])+' Mpc</td>\
                          <td >'+str(params['MJD'])+'</td>\
                          <td >'+str(d.strftime('%H:%M:%S \t %b %d, %Y UTC'))+'</td>\
                        </tr>'
        except:
            indextable += '<tr onmouseover="ChangeColor(this, true);" onmouseout="ChangeColor(this, false);" \
                        onclick="DoNav(\'Triggers/'+trig+'/'+trig+'_trigger.html\');">\
                          <td></td>\
                          <td >'+trig+'</td>\
                          <td >NA</td>\
                          <td >NA</td>\
                          <td >'+str(params['FAR'])+'</td>\
                          <td >NA</td>\
                          <td >'+str(params['MJD'])+'</td>\
                          <td >'+str(d.strftime('%H:%M:%S \t %b %d, %Y UTC'))+'</td>\
                        </tr>'
    indextable += '</tbody></table> <div class="copyright-bar bg-white"></div>'

    html = '<!DOCTYPE html> \
                <html lang="en" style="height:100%;">\
                    <head> \
                        <meta charset="utf-8"> \
                        <title>DES GW</title>\
                        <meta name="viewport" content="width=device-width, initial-scale=1.0"> \
                        <meta name="keywords" content="pinegrow, blocks, bootstrap" />\
                        <meta name="description" content="My new website" />\
                        <!-- Core CSS -->\
                        <link href="/desgw/assets/bootstrap/css/bootstrap.min.css" rel="stylesheet"> \
                        <link href="/desgw/assets/css/font-awesome.min.css" rel="stylesheet">\
                        <link href="http://fonts.googleapis.com/css?family=Open+Sans:300italic,400italic,600italic,700italic,400,300,600,700" rel="stylesheet">\
                        <link href="http://fonts.googleapis.com/css?family=Lato:300,400,700,300italic,400italic,700italic" rel="stylesheet">\
                        <!-- Style Library -->\
                        <link href="/desgw/assets/css/style-library-1.css" rel="stylesheet">\
                        <link href="/desgw/assets/css/plugins.css" rel="stylesheet">\
                        <link href="/desgw/assets/css/blocks.css" rel="stylesheet">\
                        <link href="/desgw/assets/css/custom.css" rel="stylesheet">\
                        <!--[if lt IE 9]>\
                      <script src="/desgw/assets/js/html5shiv.js"></script>\
                      <script src="/desgw/assets/js/respond.min.js"></script>\
                    <![endif]-->\
                        <script type="text/javascript">\
                    function ChangeColor(tableRow, highLight){if (highLight){\
                    tableRow.style.backgroundColor = "#7FFFD4";}\
                    else{tableRow.style.backgroundColor = "white";}}\
                  </script>\
                  <script type="text/javascript">\
                  function DoNav(theUrl){document.location.href = theUrl;}\
                  </script>\
                    </head>\
                    <body data-spy="scroll" data-target="nav">\
                        <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7.1/jquery.min.js"></script>\
                        <script type="text/javascript" src="/desgw/assets/js/smoothscroll.js"></script>\
                        <script type="text/javascript" src="/desgw/assets/js/breedjs.min.js"></script>\
                        <script type="text/javascript" src="/desgw/assets/js/jquery.sortElements.js"></script>\
                        <section id="promo-3" class="content-block promo-3 min-height-600px bg-deepocean">\
                            <div class="container text-center">\
                                <h1>LIGO GW Triggers</h1>\
                                <h2>DESGW EM Followup</h2>\
                                <a href="real-triggers.html" class="btn btn-outline btn-outline-xl outline-light">REAL Triggers</a>\
                                <a href="mock-triggers.html" class="btn btn-outline btn-outline-xl outline-light">Mock Triggers</a>\
<p></p><a href="imalive.html" class="btn btn-outline btn-outline-xl outline-light">Im Alive?</a>\
                            </div>\
                            <!-- /.container -->\
                        </section>'

    html += indextable
    html += '<script type="text/javascript" src="/desgw/assets/js/jquery-1.11.1.min.js"></script>'
    html += '<script type="text/javascript" src="/desgw/assets/js/bootstrap.min.js"></script>'
    html += '<script type="text/javascript" src="/desgw/assets/js/plugins.js"></script>'
    html += '<script src="https://maps.google.com/maps/api/js?sensor=true"></script>'
    html += '<script type="text/javascript" src="/desgw/assets/js/bskit-scripts.js"></script>\
                    </body></html>'

    if real_or_sim == 'real':

        a = open(os.path.join(webpage_dir, 'index.html'), 'w')
        a.write(html)
        a.close()
        a = open(os.path.join(webpage_dir, 'real-triggers.html'), 'w')
        a.write(html)
        a.close()
    elif real_or_sim == 'sim':
        a = open(os.path.join(webpage_dir, 'index.html'), 'w')
        a.write(html)
        a.close()
        a = open(os.path.join(webpage_dir, 'mock-triggers.html'), 'w')
        a.write(html)
        a.close()



if __name__ == '__main__':
    pass
