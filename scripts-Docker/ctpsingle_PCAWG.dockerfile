FROM smalikic/dreamchallenge_base_image
RUN git clone https://smalikic@bitbucket.org/smalikic/docker_ctpsingle_pcawg.git
RUN cp -rf docker_ctpsingle_pcawg/* /opt/CTPsingle/
