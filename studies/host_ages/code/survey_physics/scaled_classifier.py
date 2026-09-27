#!/usr/bin/env python3
"""Apply the unchanged classifier to a separate scaled-campaign workspace."""
import common

common.WORK = common.WORK / "scaled-bbc"
import classify

if __name__ == "__main__":
    classify.main()
