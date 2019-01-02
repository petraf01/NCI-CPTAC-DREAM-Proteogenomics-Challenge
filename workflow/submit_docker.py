
import synapseclient
syn = synapseclient.login()
#use docker branch of synapseclient

temp = syn.getChildren("syn8228304",includeTypes=['dockerrepo'])
for i in temp:
	syn.submit(9614155, i['id'])