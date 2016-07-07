#!/usr/bin/python

import urllib2
import urllib
import time
import sys


#run program as:
#scrapeRegDBscore.py <name of VCF file to get scores for>

def main():

    inFileName = sys.argv[1]
    rsidCol = int(sys.argv[2])

    # Iterate over entries in the input file
    with open (inFileName, 'r') as infile:  #when you use "with open" you don't have to close the file later
            for line in infile:
                if line.startswith("#"):   
                    next(infile) #skip over info and header lines
                elif line.startswith("chr"):  #the last couple lines of VCF with variant count and run time will be ignored
                    coords = line.strip("\n").split("\t")
                    # Get the regulomeDB entry based on rsid
                    pageData = getRegulomeDBDataWithText(coords[rsidCol]).strip("\n").split("\t")
                    # Append regulomeDB score to the input line and print to STDOUT
                    if len(pageData) >= 5:
                        outStr = '\t'.join(coords + [pageData[4]])
                        print(outStr)
                    else:
                        outStr = '\t'.join(coords + ['.'])
                        print(outStr)

                    
# HttpBot
# -------
# This class opens a handler and can GET and POST to a url.
#
class HttpBot:
    """an HttpBot represents one browser session, with cookies."""
    def __init__(self):
        cookie_handler= urllib2.HTTPCookieProcessor()
        redirect_handler= urllib2.HTTPRedirectHandler()
        http_handler= urllib2.HTTPHandler()
        error_handler=urllib2.HTTPErrorProcessor()
        self._opener = urllib2.build_opener(cookie_handler,redirect_handler, http_handler,error_handler)
        #urllib2.install_opener(self._opener)

    def GET(self, url):
        return self._opener.open(url).read()

    def POST(self, url, parameters):     
        user_agent = 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'
        headers = { 'User-Agent' : user_agent }
        return self._opener.open(url, urllib.urlencode(parameters))
    
# getRegulomeDBData
# -----------------
# This function takes a string of 'chr_:pos-pos\nchr_:pos-pos...' and returns a string with
# regulatory information for each site.
#
def getRegulomeDBDataWithText(coordText):
    #coords = 'chrX:55041618-55041619'
    #print (coordText)
    #print("next")
    if coordText == "": return ""
    bot = HttpBot()
    # submit coordinates to regulome db
    trials=0
    # not sure why but sometimes will get error and resubmitting fixes, so have a counter to try again and again
    while trials<20:
        try:
            x = bot.POST('http://regulome.stanford.edu/results', {'data': coordText})
            # get regulome db sid--this seems to be an id associated with the output results file
            # its needed to directly get results in text format
            pageData = x.read();
            pageData = pageData.split("<input type=\"hidden\" name=\"sid\" value=\"")
            sid = pageData[1].split("\" />")[0]
            # get output results file from regulome db using the sid
            x = bot.POST("http://regulome.stanford.edu/download", {'format':'full', 'sid':sid})  #{'format':'bed', 'sid':sid})
            #The above format can be "bed", "gff", or "full"
            pageData = x.read()
            # remove header
            pageData = pageData.split("\n", 1)
            pageData = pageData[1]
            return pageData
        except:
            # for whatever reason, errors are common so I iterate until a particular set of positions is successful or 
            # I've reached the max number of iterations
            time.sleep(0.25)
        trials+=1
    errStr = ''.join(coordText + ": Maximum attempts reached without success.\n")
    sys.stderr.write(errStr)
    return ""

#print getRegulomeDBData('chrX:55041618-55041619')

if __name__ == '__main__':
    main()
