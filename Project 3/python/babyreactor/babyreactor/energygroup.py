"""Class to contain the properties of a given energy group

"""


class EnergyGroup(object):
    """ Class to contain the properties of a given group

    """
    def __init__(self, groupnumber):
        self._groupnumber = groupnumber

    def __lt__(self, other):
        """ Returns true if neutrons from this group could eventually
            scatter into other

        """
        if self.groupnumber < other.groupnumber:
            return True
        else:
            return False

    def __gt__(self, other):
        """ Returns true if neutrons from other could eventually
            scatter into this group

        """
        if other.groupnumber > self.groupnumber:
            return True
        else:
            return False

    @property
    def groupnumber(self):
        """ The group number (index)
        """
        return self._groupnumber

    @groupnumber.setter
    def groupnumber(self, number):
        self._groupnumber = number
