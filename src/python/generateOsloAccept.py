#!/home/kmsmith/bin/python2.5
import sys
sys.path.append("/home/kmsmith/lib")
import quiet_sql
import quiet_globals
db='quiet_quality'

if (len(sys.argv)==1):
    query=open("canonical_cuts.query").read().strip()
else:
    query="SELECT q19_scan.run_id, q19_scan.run_subid, q19_timestream.module, q19_timestream.diode FROM q19_scan, q19_housekeeping, q19_timestream, q19_typeb, q19_weather WHERE %s" % getattr(quiet_globals,sys.argv[1])
data=quiet_sql.query_db(query,db)
data=set(data)
IDs = set([(row[0],row[1]) for row in data])
IDs = sorted(IDs)
nmodules=17
diodes = ["Q1","U1","U2","Q2"]
ndiodes = len(diodes)
for (run,seg) in IDs:
    accept=[]
    for module in xrange(nmodules):
        for i,diode in enumerate(diodes):
            if (run,seg,module,diode) in data:
                accept.append((module,i))
    exclusions = []
    for module in xrange(nmodules):
        for i in xrange(ndiodes):
            if (module,i) not in accept:
                exclusions.append((module,i))

    nex=len(exclusions)
    sys.stdout.write("%d %d %d " % (run,seg+1,nex))
    for ex in exclusions:
        sys.stdout.write("%d %d " % ex)
    sys.stdout.write("\n")

                
