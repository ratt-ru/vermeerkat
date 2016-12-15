#!/bin/env bash

PYTHONPATH='.' luigi --module vermeertasks --logging-conf-file logging.conf WscleanTask