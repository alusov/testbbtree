all dep clean indent tests::
	cd watchackley1 && $(MAKE) $@ && cd ..
	
doc: indent doxy

clean::
	rm -rf *~ PI* core bin/* obj/* tmp *.log
