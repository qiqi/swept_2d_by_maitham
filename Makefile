default:	bin/streamer

bin/streamer:	bin
	cd src; make ../bin/streamer

bin:
	mkdir bin
