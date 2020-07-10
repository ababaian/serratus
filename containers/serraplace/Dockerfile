FROM ubuntu:18.04
# install all the dependencies
RUN 	apt-get update &&\
	apt-get -y install wget git cmake build-essential zlib1g-dev gzip unzip flex bison &&\
	wget http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz&& tar xzvf hmmer-3.3.tar.gz && cd hmmer-3.3/ &&\
	./configure && make -j4 && make install && cd easel && make install && cd ../../ &&\
	git clone --recursive https://github.com/Pbdas/epa-ng.git &&\
	cd epa-ng && git checkout tags/v0.3.7 && git submodule update --init --recursive && make -j4 && cp bin/epa-ng /usr/local/bin && cd - &&\
	git clone https://github.com/lh3/seqtk.git && cd seqtk && make && cp seqtk /usr/local/bin && cd - &&\
	git clone --recursive https://github.com/lczech/gappa.git &&\
	cd gappa && git checkout f05e9fe && git submodule update --init --recursive && make -j4 && cp bin/gappa /usr/local/bin && cd - &&\
	cd gappa/libs/genesis/apps && git clone https://github.com/Pbdas/genesis-apps.git && make -j4 -C .. && cp ../bin/apps/genesis-apps/msa-merge ../bin/libgenesis.so /usr/local/bin && cd - &&\
	mkdir -p /serratus-data/serraplace && cd /serratus-data/serraplace && wget https://serratus-public.s3.amazonaws.com/pb/serraplace/reference/refpack.tar.gz && tar xzvf refpack.tar.gz && cd -

WORKDIR /home/serratus
COPY place.sh /home/serratus/

ENTRYPOINT ["/home/serratus/place.sh", "-d"]
