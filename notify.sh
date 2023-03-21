#!/bin/bash

TOKEN="1504507754:AAFDqHlSuYW4JEBul_8667W0uB5U6yG0Fsc"
CHAT_ID="241504810"
URL="https://api.telegram.org/bot${TOKEN}/sendMessage"

curl -s -X POST ${URL} -d chat_id=${CHAT_ID} -d text="${1}" 1>/dev/null