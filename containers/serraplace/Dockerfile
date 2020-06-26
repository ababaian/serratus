FROM ubuntu:18.04
# install all the dependencies
RUN apt-get update &&\
	apt-get -y install wget git cmake build-essential flex bison emboss hmmer &&\
	git clone --recursive https://github.com/Pbdas/epa-ng.git &&\
	cd epa-ng && git checkout tags/v0.3.7 && make -j4 && cp bin/epa-ng /usr/local/bin && cd - &&\
	git clone https://github.com/lh3/seqtk.git &&\
	cd seqtk && make && cp seqtk /usr/local/bin && cd - &&\
	wget https://github.com/amkozlov/raxml-ng/releases/download/0.9.0/raxml-ng_v0.9.0_linux_x86_64.zip &&\
	unzip raxml-ng_v0.9.0_linux_x86_64.zip raxml-ng && mv raxml-ng /usr/local/bin && rm raxml-ng_v0.9.0_linux_x86_64.zip &&\
	git clone --recursive https://github.com/lczech/gappa.git &&\
	cd gappa && git checkout dev && make -j4 && cp bin/gappa /usr/local/bin && cd - &&\
	cd gappa/libs/genesis/apps && clone https://github.com/Pbdas/genesis-apps.git && make -j4 -C .. && cp ../bin/apps/genesis-apps/msa-merge /usr/local/bin && cd -


WORKDIR /home/serratus
COPY serraplace/*sh ./


