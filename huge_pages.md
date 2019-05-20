

## Check huge pages support

Depending on the processor, there are at least two different huge page sizes on the x86_64 architecture: 2MB and 1GB. If the CPU supports 2 MB pages, it has the PSE cpuinfo flag, for 1 GB it has the PDPE1GB flag. `/proc/cpuinfo` shows whether the two flags are set. 

To check 2 MB huge page support:

    $ grep pse /proc/cpuinfo | uniq
    flags           : [...] pse [...]

To check 1 GB huge page support:

    $ grep pdpe1gb /proc/cpuinfo | uniq
    flags           : [...] pdpe1gb [...]

Reference: https://wiki.debian.org/Hugepages


## Huge pages activation

1.  First, install libhugetlbfs to ease huge pages management.
    On CentOS:

		$ sudo yum install libhugetlbfs-utils libhugetlbfs

    On Ubuntu:

        $ sudo apt-get install libhugetlbfs-dev

	`libhugetlbfs` is the library which provides easy access to hugepages of memory to any application.
	`libhugetlbfs-utils` is a set of user-space tools for configure and manage hugepage environment.
	
    Reference:
	https://paolozaino.wordpress.com/2016/10/02/how-to-force-any-linux-application-to-use-hugepages-without-modifying-the-source-code/


2.  Huge page allocation at boot time

    Allocating huge pages at boot time is less prone to failure than at run time because it is easier to find contiguous memory blocks.
    To allocate huge pages at boot time, include the huge page parameters at the `/etc/default/grub` file.
    For instance, to allocate 16 1-GiB huge pages, append the following to to the `GRUB_CMDLINE_LINUX` line:

        hugepagesz=1G hugepages=32 default_hugepagesz=1G

    Changes are applied by executing the following command: 

        $ sudo grub2-mkconfig -o /boot/grub2/grub.cfg
        Generating grub configuration file ...
        Found linux image: /boot/vmlinuz-3.10.0-693.el7.x86_64
        Found initrd image: /boot/initramfs-3.10.0-693.el7.x86_64.img
        Found linux image: /boot/vmlinuz-0-rescue-40002d05bceb4c178131ba731daaf0c4
        Found initrd image: /boot/initramfs-0-rescue-40002d05bceb4c178131ba731daaf0c4.img
        done

    And rebooting the system:

        $ sudo reboot


    To verify the parameters used at boot time:

        $ cat /proc/cmdline
        BOOT_IMAGE=/vmlinuz-3.10.0-693.el7.x86_64 root=/dev/mapper/centos-root ro 
           crashkernel=auto rd.lvm.lv=centos/root rd.lvm.lv=centos/swap rhgb quiet
           hugepagesz=1G hugepages=32 default_hugepagesz=1G
        
    To verify the huge pages allocation:
        
        $ hugeadm --pool-list
        Size  Minimum  Current  Maximum  Default
        1073741824       32       32       32        *

    In case a multiprocessor system is used (motherboard with several chips),
    we have to consider that huge pages are distributed among the processors.
    Hence, to allocate 16 1-GiB huge pages in a system composed by two processors (DP),
    we have to allocate 32 1-GB huge pages:
    
        hugepagesz=1G hugepages=32 default_hugepagesz=1G

    References:
        https://wiki.centos.org/HowTos/Grub2


3.  Huge page allocation at run time

    Allocating huge pages at runtime should be done as close to boot as possible. Memory fragments the longer the system is running and defragmenting efforts may not be able to secure enough continuous memory to allocate huge pages.

    Runtime allocations can be done in a number of ways. For instance:

		# echo 8192 > /proc/sys/vm/nr_hugepages
        $ echo 8192 | sudo tee /proc/sys/vm/nr_hugepages

    To verify the effect of the previous command:

        $ cat /sys/devices/system/node/node0/hugepages/hugepages-2048kB/nr_hugepages 
        8192

    Or:

        $ hugeadm --pool-list
              Size  Minimum  Current  Maximum  Default
           2097152     8192     8192     8192        *
        1073741824        0        0        0


    One more example, to allocate 16 1-GiB huge pages in the node 0:

        $ echo 16 | sudo tee /sys/devices/system/node/node0/hugepages/hugepages-1048576kB/nr_hugepages

    References:
	https://gist.github.com/sjenning/b6bed5bf029c9fd6f078f76b37f0a73f
    https://wiki.debian.org/Hugepages



4.  Enabling and disabling transparent huge pages (THP)

    At run time:

        $ echo always | sudo tee /sys/kernel/mm/transparent_hugepage/enabled
        $ echo always | sudo tee /sys/kernel/mm/transparent_hugepage/defrag

        $ echo never | sudo tee /sys/kernel/mm/transparent_hugepage/enabled
        $ echo never | sudo tee /sys/kernel/mm/transparent_hugepage/defrag

    To disable THP support during boot time, create the following file:

        $ sudo vi /etc/systemd/system/disable-thp.service
        
        And paste the following content:

        [Unit]
        Description=Disable Transparent Huge Pages (THP)

        [Service]
        Type=simple
        ExecStart=/bin/sh -c "echo 'never' > /sys/kernel/mm/transparent_hugepage/enabled && echo 'never' > /sys/kernel/mm/transparent_hugepage/defrag"

        [Install]
        WantedBy=multi-user.target
    
    Save the file and reload SystemD daemon:

        $ sudo systemctl daemon-reload

    Then you can start the script and enable it on boot level:

        $ sudo systemctl start disable-thp
        $ sudo systemctl enable disable-thp
        Created symlink from /etc/systemd/system/multi-user.target.wants/disable-thp.service to /etc/systemd/system/disable-thp.service.
    
    Reference:
    https://blacksaildivision.com/how-to-disable-transparent-huge-pages-on-centos

