DockerScript-dreamchallenge_base_image.txt was used to build dreamchallenge_base_image (which can be used for PCAWG as well)

DockerScript-ctpsingle_dreamchallenge.txt was used to build ctpsingle_dreamchallenge image (which was used for DREAM CHALLENGE competition)


DockerScript-ctpsingle_PCAWG.txt -- this docker needs to be built.
||
\/
--- This docker was built using the following command: 
	docker build -t ctpsingle_pcawg -f ./ctpsingle_PCAWG.dockerfile

--- It was then uploaded to docker repository following the instructions given at
	https://docs.docker.com/v17.12/docker-cloud/builds/push-images/

--- More details about installing docker can be found at:
	https://github.com/NaturalHistoryMuseum/scratchpads2/wiki/Install-Docker-and-Docker-Compose-(Centos-7)
