import sys
import pandas as pd
from time import time

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2


def _make_request(url, data):
    request = urllib2.Request(url)

    try:
        url_file = urllib2.urlopen(request)
    except urllib2.HTTPError as e:
        if e.code == 404:
            print("[NOTFOUND %d] %s" % (e.code, url))
        else:
            print("[ERROR %d] %s" % (e.code, url))

        return None

    return url_file.read().decode()


def _post_request(url, data):
    base = 'http://www.uniprot.org'
    full_url = "%s/%s" % (base, url)

    if isinstance(data, (list, tuple)):
        data = ",".join(data)

    return _make_request(full_url, data.encode())


def split_loc(s):
    try:
        return [a.split(';') for a in s.split('.') if a.split(';') != ['']]
    except AttributeError:
        return [[]]


def _set_b4_flag(s):
    return len([a for x in s['Subcellular location [CC]'] for a in x if
                   ('Cell membrane' in a or 'Secreted' in a) and ('ECO:0000269' in a or 'ECO:0000305' in a)])


def _set_b6_flag(s):
    return len([a for x in s['Subcellular location [CC]'] for a in x if
                   ('Cell membrane' in a or 'Secreted' in a) and not ('ECO:0000269' in a or 'ECO:0000305' in a)])

start = time()

base = 'http://www.uniprot.org'
# url = "uniprot/?format=xml&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&offset=50&columns=id%2centry+name%2ccomment(SUBCELLULAR+LOCATION)#"
url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
data = ['P42336', 'P60484']

location = _post_request(url, data)
location = [x.split('\t') for x in location.split('\n')]
df = pd.DataFrame(location[1:], columns=location[0])
df['Subcellular location [CC]'] = df['Subcellular location [CC]'].apply(split_loc)
df['b4'] = df.apply(_set_b4_flag, axis=1)
df['b6'] = df.apply(_set_b6_flag,axis=1)

df.to_csv('uniprot_location.csv')

print("Time taken: {}".format(time() -start))
