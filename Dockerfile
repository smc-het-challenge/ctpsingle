FROM smalikic/dreamchallenge_base_image
RUN git clone https://smalikic@bitbucket.org/smalikic/ctpsingle_dreamchallenge.git
RUN cp -rf ctpsingle_dreamchallenge/* /opt/CTPsingle/
