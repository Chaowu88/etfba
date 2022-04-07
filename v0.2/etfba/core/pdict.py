#!/usr/bin/env pyhton
# -*- coding: UTF-8 -*-


__author__ = 'Chao Wu'




class PrettyDict(dict):
	
	def __repr__(self):
		
		itemsStr = ['%s: %s' % (key, value) for key, value in sorted(self.items())]
		
		return '\n'.join(itemsStr)
		
		
		

class SmartDict(PrettyDict):
	
	def __init__(self, *args, **kwargs):
		
		super().__init__(*args, **kwargs)
		self.caller = None
		
	
	def __getitem__(self, key):
		
		item = super().__getitem__(key)
		if self.caller != None:
			item.host = self.caller
			
		return item	
		
		
		