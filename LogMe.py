#!/usr/bin/env python

import time

class LogMe:
	def __init__(self,logFileName):
		self.logFileName=logFileName
		self.data=[]

	def __del__(self):
		self.SaveLog(self.logFileName)


	def SaveLog(self,logFileName):
		with open(logFileName,'w') as outfile:
			outfile.write(time.strftime('%c')+'\n')
			for cline in self.data:
				outfile.write(cline+'\n')


	def log(self,*args):
		cstring=''
		cstring=' '.join(str(obj) for obj in args)
		self.data.append(cstring)

def main():
	log=LogMe()
	log.log("This is a test")
	log.log("and this","is a test",2)
	log.log("The end")

if __name__ == "__main__":
    main()