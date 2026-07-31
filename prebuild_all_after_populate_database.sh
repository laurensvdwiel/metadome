  #!/bin/bash
  set -e
  COMPOSE_FILE="docker-compose.yml"

  # Optional: pull a custom compose file off the front, forward everything else to the module.
  if [[ "$1" == "-c" || "$1" == "--compose-file" ]]; then
      COMPOSE_FILE="$2"; shift 2
  fi

  # Reuse the (possibly live) stack: ensure only db is up. Never `down`, never `-v`.
docker compose -f "$COMPOSE_FILE" up -d --wait db

  ( while true; do
      docker stats --no-stream --format 'table {{.Name}}\t{{.CPUPerc}}\t{{.MemUsage}}\t{{.BlockIO}}' >> prebuild_stats.log 2>/dev/null
      sleep 15
    done ) &
  trap 'kill $! 2>/dev/null' EXIT

  # --no-deps: synchronous prebuild needs no broker/workers/mailserver.
  # --rm: removes only this one-shot container on exit (foreground, so Ctrl+C reaches it).
  docker compose -f "$COMPOSE_FILE" run --rm --no-deps \
    app python -m metadome.prebuild_all "$@"