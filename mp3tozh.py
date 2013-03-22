#encoding=utf8
from path import path
import subprocess, time
cwd = path.getcwd()
for mp3 in cwd.walkfiles('*.mp3'):
	subprocess.call(['mid3iconv', '-e', 'GBK', mp3])
