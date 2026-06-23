BEGIN {
    FS=", "
    basilisk = ENVIRON["BASILISK"]
}

function valid_ref(s) {
    return s ~ /^[A-Za-z0-9_.:-]+$/
}

function valid_basilisk(s) {
    return s ~ /^[A-Za-z0-9_\/.+:-]+$/
}

/@hal{.*,.*}/ {
    key = $1
    hal_id = $2
    gsub(/}/, "", hal_id)
    gsub(/.*{/, "", key)

    if (!valid_ref(key) || !valid_ref(hal_id)) {
        print "@misc{error, title={invalid HAL reference}}"
        next
    }

    if (!valid_basilisk(basilisk)) {
        print "@misc{error, title={invalid Darcsit configuration}}"
        next
    }

    cmd = "bash " basilisk "/darcsit/gethal.sh " hal_id
    while ((cmd | getline line) > 0) {
        if (line ~ /^[[:space:]]*@[^{]+{[^,]*,/) {
            sub(/{[^,]*,/, "{" key ",", line)
        }
        print line
    }
    close(cmd)
    next
}

{ print $0 }
