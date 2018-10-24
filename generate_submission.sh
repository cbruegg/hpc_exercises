#!/bin/bash

zip submission.zip $(git ls-tree -r HEAD --name-only | awk 1 ORS=' ')