if ! grep -q "cd /home/vagrant/build" ~/.bashrc ; then
    echo "cd /home/vagrant/build" >> ~/.bashrc
fi
sudo apt-get update && \
sudo apt-get install debootstrap
