import multiprocessing

workers = 20
bind = 'unix:flaskrest.sock'
umask = 0o007
reload = True

#logging
accesslog = '-'
errorlog = '-'
