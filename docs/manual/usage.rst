Usage
=====

Typical usage is of course heavily tied to data and the `E200 <http://e200.readthedocs.org/>`_ package.

First, open up some data and select an image::

    import E200
    filename = 'nas/nas-li20-pm00/E217/2014/20140601/E217_13123/E217_13123.mat'
    data     = E200.E200_load_data(filename)
    cam      = data.rdrill.data.raw.images.ELANEX
    shot_ind = 39
    uids     = cam.UID[shot_ind]
    imgs     = E200.E200_load_images(cam, UID=uids)

Then, feed the image into things
