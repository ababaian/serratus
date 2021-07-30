##### Builder Image #######
FROM base_img_serratus:latest
# install all the dependencies
RUN 	apk add cmake bison flex zlib-dev && git clone --recursive https://github.com/Pbdas/epa-ng.git &&\
	cd epa-ng && git checkout tags/v0.3.7 && git submodule update --init --recursive && make -j4 && cd - &&\
	git clone https://github.com/lh3/seqtk.git && cd seqtk && make && cd - &&\
	git clone --recursive https://github.com/lczech/gappa.git &&\
	cd gappa && git checkout f05e9fe && git submodule update --init --recursive && make -j4 && cd - &&\
	cd gappa/libs/genesis/apps && git clone https://github.com/Pbdas/genesis-apps.git && make -j4 -C .. && cd -

#######  Final Image #########

FROM base_img_serratus:latest
WORKDIR /home/serratus
COPY --from=builder epa-ng/bin/epa-ng  /usr/local/bin
COPY --from=builder seqtk/seqtk /usr/local/bin
COPY --from=builder gappa/libs/genesis/bin/gappa gappa/libs/genesis/bin/libgenesis.so /usr/local/bin
COPY place.sh /home/serratus/

ENTRYPOINT ["/home/serratus/place.sh", "-d"]
