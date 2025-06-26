#!/usr/bin/env bash

DOCKER_IMAGE=quay.io/bgruening/galaxy
PORT=8080
tool_dir="$(dirname $(dirname $(readlink -f "$0")))/tools"

status() {
    if [[ $(sudo lsof -i :$PORT) ]]; then
        echo "Galaxy is running !"
    else
        echo "Galaxy is not running !"
    fi
}

start() {


        # launching docker compose in detached mode
        echo "Launching Galaxy in detached mode."
        docker run \
            -d \
            -p 8080:80 \
            -p 8021:21 \
            -p 8022:22 \
            -v $tool_dir:/local_tools \
            -e GALAXY_CONFIG_TOOL_CONFIG_FILE=/etc/galaxy/tool_conf.xml,/local_tools/tool_conf.xml \
            $DOCKER_IMAGE
        echo "Galaxy started !"
}

stop() {
        docker stop $(docker ps | grep galaxy | awk '{print $1}')
        echo "Galaxy stopped !"
}

case "$1" in
    'start')
            start
            ;;
    'stop')
            stop
            ;;
    'restart')
            stop
            start
            ;;
    'status')
            status
            ;;
    *)
            echo
            echo "Usage: $0 { start | stop | restart | status }"
            echo
            exit 1
            ;;
esac

exit 0

